import Base: size
import LinearSolver: initkaczmarz, dot_with_matrix_row, kaczmarz_update!

export reconstructionSeparate, reconstructionFFJoint, FFOperator

# this wrapper should maybe be removed
function reconstruction{T<:MPIFile}(bSF::MPIFile, bMeas::Vector{T}, freq; kargs...)
  reconstructionFFJoint(bSF, bMeas, freq; kargs...)
end

# Necessary for Multi-System-Matrix FF reconstruction
voxelSize(bSF::Vector{T}) where T<:MPIFile = voxelSize(bSF[1])
sfGradient(bSF::Vector{T},dim) where T<:MPIFile = sfGradient(bSF[1],dim)
generateHeaderDict(bSF::Vector{T},bMeas::MPIFile) where T<:MPIFile =
   generateHeaderDict(bSF[1],bMeas)

function reconstructionFFJoint(bSF, bMeas::MPIFile, freq;
            frames=nothing, bEmpty=nothing, OverscanSF=[0,0,0],OffsetFF=[0.002,0.002,0.002],
            nAverages=1,denoiseWeight=0, loadas32bit=false, FFPos = ffPos(bMeas),
            spectralLeakageCorrection=true, kargs...)

  consistenceCheck(bSF, bMeas)

  bgcorrection = (bEmpty != nothing)


  periodsSortedbyFFPos = unflattenOffsetFieldShift(FFPos)

  FFPos = FFPos[:,periodsSortedbyFFPos[:,1]]

  FFOp = FFOperator(bSF,bMeas,freq,bgcorrection,OverscanSF=OverscanSF,
                    OffsetFF=OffsetFF,denoiseWeight=denoiseWeight,FFPos=FFPos)

  L = numScans(bMeas)
  (frames==nothing) && (frames=collect(1:L))
  nFrames=length(frames)

  uTotal = getMeasurementsFD(bMeas,frequencies=freq, frames=frames, numAverages=nAverages,
                             spectralLeakageCorrection=spectralLeakageCorrection)

  uTotal=uTotal[:,periodsSortedbyFFPos,:]
  uTotal=mean(uTotal,3)

  # Here we call a regular reconstruction function
  c = reconstruction(FFOp,uTotal,(shape(FFOp.grid)...,); kargs...)

  pixspacing = voxelSize(bSF) ./ sfGradient(bMeas,3) .* sfGradient(bSF,3)
  offset = fieldOfViewCenter(FFOp.grid)  .- 0.5.*fieldOfView(FFOp.grid) .+ 0.5.*spacing(FFOp.grid)

  im = AxisArray(c, (:x,:y,:z,:time), tuple(pixspacing...,dfcycle(bMeas)),
                                      tuple(offset...,0.0))

  imMeta = ImageMeta(im,generateHeaderDict(bSF,bMeas))
  return imMeta
end


# FFOperator is a type that acts as the MPI system matrix but exploits
# its sparse structure.
# Its very important to keep this type typestable
type FFOperator{V<:AbstractMatrix, T<:Positions}
  S::Vector{V}
  grid::T
  N::Int
  M::Int
  RowToPatch::Vector{Int}
  xcc::Vector{Vector{Int}}
  xss::Vector{Vector{Int}}
  nPatches::Int
  patchToSMIdx::Vector{Int}
end


function FFOperator(SF::MPIFile, bMeas, freq::Vector{Int64}, bgcorrection::Bool; kargs...)
  return FFOperator([SF], bMeas, freq, bgcorrection; kargs...)
end

function findNearestPatch(ffPosSF, FFPos)
  idx = -1
  minDist = 1e20
  for (l,FFPSF) in enumerate(ffPosSF)
    dist = norm(FFPSF.-FFPos)
    if dist < minDist
      minDist = dist
      idx = l
    end
  end
  if idx < 0
    error("Something went wrong")
  end
  return idx
end

function FFOperator(SFs::Vector, bMeas, freq::Vector{Int64}, bgcorrection::Bool;
                    denoiseWeight=0, FFPos=zeros(0,0), patchMirroring = false, kargs...)

  println("Load SF")
  # Maybe introduce interpolation here ( if grids do not fit )
  S = [getSF(SF,freq,nothing,"kaczmarz", bgcorrection=bgcorrection) for SF in SFs]
  numPatches = size(FFPos,2)
  M = size(S[1],1)
  RowToPatch = kron(collect(1:numPatches), ones(Int,M))

  if !patchMirroring
    ffPosSF = [vec(ffPos(SF)) for SF in SFs]
    ffPosSFAbs = [vec(abs.(ffPos(SF))) for SF in SFs]
    
    positions = CartesianGridPositions[]
    patchToSMIdx = zeros(Int,numPatches) 
    
    for k=1:numPatches
      idx = findNearestPatch(ffPosSF, FFPos[:,k])
      SF = SFs[idx]
      diffFFPos = ffPosSF[idx] .- FFPos[:,k]
          
      push!(positions, CartesianGridPositions(calibSize(SF),calibFov(SF),calibFovCenter(SF).-diffFFPos))
      patchToSMIdx[k] = idx
    end 
  else # patchMirroring
    ffPosSF = [vec(ffPos(SF)) for SF in SFs]
    ffPosSFAbs = [vec(abs.(ffPos(SF))) for SF in SFs]
    
    positions = CartesianGridPositions[]
    patchToSMIdx = zeros(Int,numPatches) 
    
    for k=1:numPatches
      
      idx = findfirst(x -> isapprox(x,FFPos[:,k]),ffPosSF)
      if idx == 0
        idx = findfirst(x -> isapprox(x,abs.(FFPos[:,k])),ffPosSFAbs)
      end
      SF = SFs[idx]
      if idx > 0
        signs = [isapprox(ffPosSF[idx][d],FFPos[d,k]) ? 1 : -1 for d=1:3]
        diffFFPos = ffPosSF[idx] .- FFPos[:,k]
          
        push!(positions, CartesianGridPositions(calibSize(SF),calibFov(SF),calibFovCenter(SF).-diffFFPos, abs.(signs)))
        patchToSMIdx[k] = idx
      else
        error("Did not find a suitable Calibration Scan!  $(FFPos[:,k]) \n $(ffPosSFAbs)")
      end
    end  
  end

  println("Calc Grids")
  recoGrid = CartesianGridPositions(positions)

  println("Calc LUT")
  xcc, xss = calculateLUT(S, patchToSMIdx, positions, recoGrid, numPatches)  
  
  println("Finished")
  return FFOperator(S, recoGrid, length(recoGrid), M*numPatches,
             RowToPatch, xcc, xss, numPatches, patchToSMIdx)
end

function calculateLUT(S, patchToSMIdx, positions, recoGrid, numPatches)
  xss = Vector{Int}[]
  xcc = Vector{Int}[]
  for k=1:numPatches
    N = size(S[patchToSMIdx[k]],2)
    push!(xss, collect(1:N))
    xc = zeros(Int64,N)
    for n=1:N
      xc[n] = posToLinIdx(recoGrid,positions[k][n])
    end
    push!(xcc, xc)
  end
  return xcc, xss
end

function size(FFOp::FFOperator,i::Int)
  if i==2
    return FFOp.N
  elseif i==1
    return FFOp.M
  else
    error("bounds error")
  end
end

length(FFOp::FFOperator) = size(FFOp,1)*size(FFOp,2)

### The following is intended to use the standard kaczmarz method ###

function calculateTraceOfNormalMatrix(Op::FFOperator, weights)
  if length(Op.S) == 1
    trace = calculateTraceOfNormalMatrix(Op.S[1],weights)
    trace *= Op.nPatches*prod(Op.PixelSizeSF)/prod(Op.PixelSizeC)
  else
    trace = sum([calculateTraceOfNormalMatrix(S,weights) for S in Op.S])
    trace *= prod(Op.PixelSizeSF)/prod(Op.PixelSizeC)
  end
  return trace
end

setlambda(::FFOperator, ::Real) = nothing

function dot_with_matrix_row{T}(Op::FFOperator, x::AbstractArray{T}, k::Integer)
  p = Op.RowToPatch[k]
  xs = Op.xss[p]
  xc = Op.xcc[p]

  j = mod1(k,div(Op.M,Op.nPatches))
  A = Op.S[Op.patchToSMIdx[p]]

  return dot_with_matrix_row_(A,x,xs,xc,j)
end

function dot_with_matrix_row_{T}(A::AbstractArray{T},x,xs,xc,j)
  tmp = zero(T)
  @simd  for i = 1:length(xs)
     @inbounds tmp += A[j,xs[i]]*x[xc[i]]
  end
  tmp
end

function kaczmarz_update!(Op::FFOperator, x::AbstractArray, k::Integer, beta)
  p = Op.RowToPatch[k]
  xs = Op.xss[p]
  xc = Op.xcc[p]

  j = mod1(k,div(Op.M,Op.nPatches))
  A = Op.S[Op.patchToSMIdx[p]]

  kaczmarz_update_!(A,x,beta,xs,xc,j)
end

function kaczmarz_update_!(A,x,beta,xs,xc,j)
  @simd for i = 1:length(xs)
    @inbounds x[xc[i]] += beta* conj(A[j,xs[i]])
  end
end

function initkaczmarz(Op::FFOperator,λ,weights::Vector)
  T = typeof(real(Op.S[1][1]))
  denom = zeros(T,Op.M)
  rowindex = zeros(Int64,Op.M)

  MSub = div(Op.M,Op.nPatches)

  if length(Op.S) == 1
    for i=1:MSub
      s² = rownorm²(Op.S[1],i)*weights[i]^2
      if s²>0
        for l=1:Op.nPatches
          k = i+MSub*(l-1)
          denom[k] = weights[i]^2/(s²+λ)
          rowindex[k] = k
        end
      end
    end
  else
    for l=1:Op.nPatches
      for i=1:MSub
        s² = rownorm²(Op.S[Op.patchToSMIdx[l]],i)*weights[i]^2
        if s²>0
          k = i+MSub*(l-1)
          denom[k] = weights[i]^2/(s²+λ)
          rowindex[k] = k
        end
      end
    end
  end

  denom, rowindex
end
