import Base: size, eltype
import RegularizedLeastSquares: initkaczmarz, dot_with_matrix_row, kaczmarz_update!

export reconstructionMultiPatch, MultiPatchOperator


# Necessary for Multi-System-Matrix FF reconstruction
voxelSize(bSF::MultiMPIFile) = voxelSize(bSF[1])
sfGradient(bSF::MultiMPIFile,dim) = sfGradient(bSF[1],dim)
generateHeaderDict(bSF::MultiMPIFile,bMeas::MPIFile) =
   generateHeaderDict(bSF[1],bMeas)


function reconstructionMultiPatch(bSF, bMeas::MPIFile;
  minFreq=0, maxFreq=1.25e6, SNRThresh=-1,maxMixingOrder=-1, numUsedFreqs=-1, sortBySNR=false, recChannels=1:numReceivers(bMeas), kargs...)
  Base.depwarn("`reconstructionMultiPatch` is deprecated. Use `reconstruct(\"MultiPatch\", bMeas; kwargs...)` or another fitting MultiPatch algorithm instead", :reconstructionMultiPatch)

  freq = filterFrequencies(bSF,minFreq=minFreq, maxFreq=maxFreq,recChannels=recChannels, SNRThresh=SNRThresh, numUsedFreqs=numUsedFreqs)
  freq = sortFrequencies(freq, bSF, numPeriodGrouping = 1, sortBySNR = sortBySNR)

  @debug "selecting $(length(freq)) frequencies"

  return reconstructionMultiPatch(bSF, bMeas, freq; kargs...)
end

function reconstructionMultiPatch(bSF, bMeas::MPIFile, freq;
            frames=nothing, bEmpty=nothing, nAverages=1, numAverages=nAverages,
	    spectralLeakageCorrection=true, voxelsize=voxelSize(bSF), tfCorrection=rxHasTransferFunction(bSF),
	    kargs...)

  bgCorrection = (bEmpty != nothing)

  FFOp = MultiPatchOperatorHighLevel(bSF, bMeas, freq, bgCorrection;
                    kargs... )

  L = acqNumFGFrames(bMeas)
  (frames==nothing) && (frames=collect(1:L))
  nFrames=length(frames)

  uTotal_ = getMeasurementsFD(bMeas,frequencies=freq, frames=frames, numAverages=numAverages,
                             spectralLeakageCorrection=spectralLeakageCorrection, tfCorrection=tfCorrection)

  periodsSortedbyFFPos = unflattenOffsetFieldShift(ffPos(bMeas))
  uTotal = similar(uTotal_,size(uTotal_,1),length(periodsSortedbyFFPos),size(uTotal_,3))

  for k=1:length(periodsSortedbyFFPos)
      uTotal[:,k,:] = mean(uTotal_[:,periodsSortedbyFFPos[k],:], dims=2)
  end

  # subtract background measurement
  if bEmpty != nothing
    uEmpty = getMeasurementsFD(bEmpty, frequencies=freq, numAverages=1, spectralLeakageCorrection=true, tfCorrection=tfCorrection)
    numFrames = acqNumPeriodsPerPatch(bEmpty)
    bgFrames = [1+(i-1)*numFrames:i*numFrames for i=1:acqNumPatches(bEmpty)]
    for i=1:acqNumPatches(bMeas)
      uTotal[:,i,:] = uTotal[:,i,:] .- mean(uEmpty[:,bgFrames[i],:],dims=(2,3))
    end
  end

  # Here we call a regular reconstruction function
  c_ = reconstruction(FFOp, uTotal; kargs...)
  c = reshape(c_, shape(FFOp.grid)..., :)

  # calculate axis
  shp = size(c)
  pixspacing = (voxelsize ./ sfGradient(bMeas,3) .* sfGradient(bSF,3)) * 1000u"mm"
  offset = (fieldOfViewCenter(FFOp.grid) .- 0.5.*fieldOfView(FFOp.grid) .+ 0.5.*spacing(FFOp.grid)) * 1000u"mm"

  # TODO does this provide the correct value in the multi-patch case?
  dtframes = acqNumAverages(bMeas)*dfCycle(bMeas)*numAverages*1u"s"
  # create image
  c = reshape(c,1,size(c)...)
  im = makeAxisArray(c, pixspacing, offset, dtframes)
  imMeta = ImageMeta(im,generateHeaderDict(bSF,bMeas))
  return imMeta
end

export AbstractMultiPatchOperatorParameter
abstract type AbstractMultiPatchOperatorParameter <: AbstractMPIRecoParameters end
# MultiPatchOperator is a type that acts as the MPI system matrix but exploits
# its sparse structure.
# Its very important to keep this type typestable
export AbstractMultiPatchOperator, MultiPatchOperator, DenseMultiPatchOperator
abstract type AbstractMultiPatchOperator{T, V} <: AbstractArray{T, 2} end
mutable struct MultiPatchOperator{T, V <: AbstractMatrix{T}, U<:Positions, I <: Integer, vecI <: AbstractVector{I}, matI <: AbstractMatrix{I}} <: AbstractMultiPatchOperator{T, V}
  S::Vector{V}
  grid::U
  N::Int64
  M::Int64
  RowToPatch::vecI
  xcc::Vector{vecI}
  xss::Vector{vecI}
  sign::matI
  nPatches::I
  patchToSMIdx::vecI
end

mutable struct DenseMultiPatchOperator{T, V <: AbstractArray{T, 3}, U<:Positions, I <: Integer, vecI <: AbstractVector{I}, matI <: AbstractMatrix{I}} <: AbstractMultiPatchOperator{T, V}
  S::V
  grid::U
  N::Int64
  M::Int64
  RowToPatch::vecI
  xcc::matI
  xss::matI
  sign::matI
  nPatches::I
  patchToSMIdx::vecI
end
Adapt.adapt_structure(::Type{arrT}, op::MultiPatchOperator{T, <:arrT}) where {T, arrT} = op
function Adapt.adapt_structure(arr, op::MultiPatchOperator)
  S = adapt.(arr, op.S)
  grid = op.grid
  N = op.N
  M = op.M
  RowToPatch = adapt(arr, op.RowToPatch)
  xcc = adapt.(arr, op.xcc)
  xss = adapt.(arr, op.xss)
  sign = adapt(arr, op.sign)
  nPatches = op.nPatches
  patchToSMIdx = adapt(arr, op.patchToSMIdx)
  return MultiPatchOperator(S, grid, N, M, RowToPatch, xcc, xss, sign, nPatches, patchToSMIdx)
end
LinearOperators.storage_type(op::MultiPatchOperator) = LinearOperators.storage_type(first(op.S))
LinearOperators.storage_type(op::DenseMultiPatchOperator) = typeof(similar(op.S, 0))

function Base.hash(op::AbstractMultiPatchOperator, h::UInt64)
  h = hash(typeof(op), h)
  for field in fieldnames(typeof(op))
    h = hash(getfield(op, field), h)
  end
  return h
end

function Base.convert(::Type{DenseMultiPatchOperator}, op::MultiPatchOperator)
  S = stack(op.S)
  xcc = stack(op.xcc)
  xss = stack(op.xss)
  return DenseMultiPatchOperator(S, op.grid, op.N, op.M, op.RowToPatch, xcc, xss, op.sign, op.nPatches, op.patchToSMIdx)
end

eltype(FFOp::MultiPatchOperator) = eltype(FFOp.S[1])
eltype(FFOp::DenseMultiPatchOperator) = eltype(FFOp.S)

function MultiPatchOperatorHighLevel(SF::MPIFile, bMeas, freq, bgCorrection::Bool; kargs...)
  return MultiPatchOperatorHighLevel(MultiMPIFile([SF]), bMeas, freq, bgCorrection; kargs...)
end

function MultiPatchOperatorHighLevel(bSF::MultiMPIFile, bMeas, freq, bgCorrection::Bool;
        FFPos = zeros(0,0), FFPosSF = zeros(0,0), tfCorrection=rxHasTransferFunction(bSF), kargs...)

  FFPos_ = ffPos(bMeas)

  periodsSortedbyFFPos = unflattenOffsetFieldShift(FFPos_)
  idxFirstPeriod = getindex.(periodsSortedbyFFPos,1)
  FFPos_ = FFPos_[:,idxFirstPeriod]

  if length(FFPos) > 0
    FFPos_[:] = FFPos
  end

  if length(FFPosSF) == 0
    L = length(ffPos(bSF[1]))
    FFPosSF_ = [vec(ffPos(SF))[l] for l=1:L, SF in bSF] #[vec(ffPos(SF)) for SF in bSF]
  else
    FFPosSF_ = FFPosSF #[vec(FFPosSF[:,l]) for l=1:size(FFPosSF,2)]
  end

  gradient = acqGradient(bMeas)[:,:,1,idxFirstPeriod]

  # Use transfer function correction only if all system matrices and measurements have one     
  if !rxHasTransferFunction(bSF) && !rxHasTransferFunction(bMeas) && tfCorrection
    @warn "One or multiple system matrices or the measurement don't have a transfer function. No transfer function correction will be applied."
    tfCorrection = false
  end
  FFOp = MultiPatchOperator(bSF, bMeas, freq, bgCorrection;
                  FFPos = FFPos_,
                  gradient = gradient,
                  FFPosSF = FFPosSF_,
		              tfCorrection = tfCorrection, kargs...)
  return FFOp
end

function process(::Type{<:AbstractMPIRecoAlgorithm}, params::AbstractMultiPatchOperatorParameter, bSF::MultiMPIFile, freq, gradient, FFPos, FFPosSF)
  @info "Loading Multi Patch operator"
  return MultiPatchOperator(bSF, freq; toKwargs(params)..., FFPos = FFPos, FFPosSF = FFPosSF, gradient = gradient)
end 
function process(::Type{<:AbstractMPIRecoAlgorithm}, params::AbstractMultiPatchOperatorParameter, op::AbstractMultiPatchOperator, arrayType::Type{<:AbstractArray})
  return adapt(arrayType, op)
end

function MultiPatchOperator(SF::MPIFile, freq, bgCorrection::Bool; kargs...)
  return MultiPatchOperator(MultiMPIFile([SF]), freq, bgCorrection; kargs...)
end

function findNearestPatch(ffPosSF, FFPos, gradientSF, gradient)
  idx = -1
  minDist = 1e20
  for l = 1:size(ffPosSF,2)
    if gradientSF[l][:,:,1,1] == gradient
      dist = norm(ffPosSF[:,l].-FFPos)
      if dist < minDist
        minDist = dist
        idx = l
      end
    end
  end
  if idx < 0
    error("Something went wrong")
  end
  return idx
end

MultiPatchOperator(SFs, meas, freq, bgCorr; kwargs...) = MultiPatchOperator(SFs, freq; kwargs..., bgCorrection = bgCorr)
function MultiPatchOperator(SFs::MultiMPIFile, freq;
        mapping=zeros(0),fov=hcat(calibFov.(SFs)...), gridsize=hcat(calibSize.(SFs)...), center = hcat(calibFovCenter.(SFs)...), kargs...)
  if length(mapping) > 0
    return MultiPatchOperatorExpliciteMapping(SFs,freq; mapping=mapping, gridsize=gridsize, fov=fov, SFGridCenter=center, kargs...)
  else
    if any(fov .> hcat(calibFov.(SFs)...))      
        mapping=collect(1:length(SFs))
        @warn "You try to performe a system matrix extrapolation on multi-patch data without giving an explicit mapping.
Thus, the mapping is automatically set to $mapping."
        return MultiPatchOperatorExpliciteMapping(SFs,freq; mapping=mapping, gridsize=gridsize, fov=fov, SFGridCenter=center, kargs...)
    else
      return MultiPatchOperatorRegular(SFs,freq; kargs...)
    end
  end
end

export ExplicitMultiPatchParameter
Base.@kwdef struct ExplicitMultiPatchParameter <: AbstractMultiPatchOperatorParameter
  bgCorrection::Bool = false
  tfCorrection::Bool = true
  SFGridCenter::AbstractArray = zeros(0,0)
  systemMatrices::Union{Nothing, AbstractArray} = nothing
  mapping::Vector{Int64}
  gridsize::Union{Nothing, AbstractArray} = nothing
  fov::Union{Nothing, AbstractArray} = nothing
  grid::Union{Nothing, RegularGridPositions} = nothing
end
function MultiPatchOperatorExpliciteMapping(SFs::MultiMPIFile, freq; bgCorrection::Bool,
                    denoiseWeight=0, FFPos=zeros(0,0), FFPosSF=zeros(0,0),
                    gradient=zeros(0,0,0),
                    roundPatches = false,
                    SFGridCenter = hcat(calibFovCenter.(SFs)...),
                    systemMatrices = nothing,
                    mapping=zeros(0),
		                gridsize = hcat(calibSize.(SFs)...),
                    fov = hcat(calibFov.(SFs)...),
                    grid = nothing,
		                tfCorrection = rxHasTransferFunction(SFs),
		                kargs...)

                    
  @debug "Loading System matrix"
  numPatches = size(FFPos,2)
  M = length(freq)
  RowToPatch = kron(collect(1:numPatches), ones(Int,M))

  gridsize = isnothing(gridsize) ? hcat(calibSize.(SFs)...) : gridsize
  fov = isnothing(fov) ? hcat(calibFov.(SFs)...) : fov
  if length(SFGridCenter) == 0
    SFGridCenter = zeros(3,length(SFs))
    for l=1:length(SFs)
      SFGridCenter[:,l] = calibFovCenter(SFs[l])
    end
  end

  if systemMatrices == nothing
      S=AbstractMatrix[]; gridS=RegularGridPositions[];
      for i=1:length(SFs)
        SF_S,SF_gridS = getSF(SFs[i],freq,nothing,"Kaczmarz", bgCorrection=bgCorrection, tfCorrection=tfCorrection,gridsize=gridsize[:,i],fov=fov[:,i],center=SFGridCenter[:,i])
        push!(S,SF_S)
        push!(gridS,SF_gridS)
      end
  else
    if grid == nothing
      gridS = [getSF(SFs[i],freq,nothing,"Kaczmarz", bgCorrection=bgCorrection,tfCorrection=tfCorrection,gridsize=gridsize[:,i],fov=fov[:,i],center=SFGridCenter[:,i])[2] for i in 1:length(SFs)]
    else
      gridS = grid
    end
    S = systemMatrices
  end

  gradientSF = [acqGradient(SF) for SF in SFs]

  grids = RegularGridPositions[]
  patchToSMIdx = mapping


  sign = ones(Int, M, numPatches)

  # We first check which system matrix fits best to each patch. Here we use only
  # those system matrices where the gradient matches. If the gradient matches, we take
  # the system matrix with the closes focus field shift
  for k=1:numPatches
    idx = mapping[k]
    SF = (length(SFs) == 1) ? SFs[1] : SFs[idx] # systemMatrices can contain more system matrices than SFs

    diffFFPos = FFPosSF[:,idx] .- FFPos[:,k]

    # use the grid parameters of the SM or use given parameter
    calibsizeTmp = gridS[idx].shape
    calibfovTmp = gridS[idx].fov

    push!(grids, RegularGridPositions(calibsizeTmp,calibfovTmp,SFGridCenter[:,idx].-diffFFPos))
  end

  #if SMextrapolation != nothing
  #  for k=1:numPatches
  #    S[k],grids[k] = extrapolateSM(S[k],grids[k],SMextrapolation)
  #  end
  #end

  # We now know all the subgrids for each patch, if the corresponding system matrix would be taken as is
  # and if a possible focus field missmatch has been taken into account (by changing the center)
  @debug "Calculate Reconstruction Grid"
  if grid == nothing
    recoGrid = RegularGridPositions(grids)
  else
    recoGrid = grid
    # Test whether given grid covers all subgrids
    recoGridTest = RegularGridPositions(grids)
    testMin = round.(recoGridTest.center,digits=3) - round.(recoGridTest.fov./2,digits=3)
    testMax = round.(recoGridTest.center,digits=3) + round.(recoGridTest.fov./2,digits=3)
    gridMin = round.(recoGrid.center,digits=3) - round.(recoGrid.fov./2,digits=3)
    gridMax = round.(recoGrid.center,digits=3) + round.(recoGrid.fov./2,digits=3)
    if false in ((gridMin .<= testMin) .* (gridMax .>= testMax))
        @warn "A larger grid is used for reconstruction, since the given grid does not cover all subgrids: $recoGrid ⊅  $recoGridTest"
        recoGrid = recoGridTest
    end
  end


  # Within the next loop we will refine our grid since we now know our reconstruction grid
  for k=1:numPatches
    #idx = mapping[k]
    #SF = SFs[idx]

    issubgrid = isSubgrid(recoGrid,grids[k])
    if !issubgrid
      grids[k] = deriveSubgrid(recoGrid, grids[k])
    end
  end
  @debug "Use $(length(S)) patches"

  @debug "Calculate LUT"
  # now that we have all grids we can calculate the indices within the recoGrid
  xcc, xss = calculateLUT(grids, recoGrid)

  return MultiPatchOperator{eltype(first(S)), reduce(promote_type, typeof.(S)), typeof(recoGrid), eltype(patchToSMIdx), typeof(patchToSMIdx), typeof(sign)}(S, recoGrid, length(recoGrid), M*numPatches,
             RowToPatch, xcc, xss, sign, numPatches, patchToSMIdx)
end

export RegularMultiPatchOperatorParameter
Base.@kwdef struct RegularMultiPatchOperatorParameter <: AbstractMultiPatchOperatorParameter
  bgCorrection::Bool = false
  denoiseWeight::Float64=0.0
  roundPatches::Bool = false
  tfCorrection::Bool = true
end
function MultiPatchOperatorRegular(SFs::MultiMPIFile, freq; bgCorrection::Bool,
                    denoiseWeight=0, gradient=zeros(0,0,0),
                    roundPatches = false, FFPos=zeros(0,0), FFPosSF=zeros(0,0),
		    tfCorrection = true,
                    kargs...)

  @debug "Loading System matrix"
  numPatches = size(FFPos,2)
  M = length(freq)
  RowToPatch = kron(collect(1:numPatches), ones(Int,M))

  S = AbstractMatrix[]
  SOrigIdx = Int[]
  SIsPlain = Bool[]

  gradientSF = [acqGradient(SF) for SF in SFs]

  grids = RegularGridPositions[]
  matchingSMIdx = zeros(Int,numPatches)
  patchToSMIdx = zeros(Int,numPatches)

  # We first check which system matrix fits best to each patch. Here we use only
  # those system matrices where the gradient matches. If the gradient matches, we take
  # the system matrix with the closes focus field shift
  for k=1:numPatches
    idx = findNearestPatch(FFPosSF, FFPos[:,k], gradientSF, gradient[:,:,k])

    SF = SFs[idx]

    if isapprox(FFPosSF[:,idx],FFPos[:,k])
      diffFFPos = zeros(3)
    else
      diffFFPos = FFPosSF[:,idx] .- FFPos[:,k]
    end
    push!(grids, RegularGridPositions(calibSize(SF),calibFov(SF),calibFovCenter(SF).-diffFFPos))
    matchingSMIdx[k] = idx
  end

  # We now know all the subgrids for each patch, if the corresponding system matrix would be taken as is
  # and if a possible focus field missmatch has been taken into account (by changing the center)
  @debug "Calculate Reconstruction Grid"
  recoGrid = RegularGridPositions(grids)

  # Within the next loop we will refine our grid since we now know our reconstruction grid
  for k=1:numPatches
    idx = matchingSMIdx[k]
    SF = SFs[idx]

    issubgrid = isSubgrid(recoGrid,grids[k])
    if !issubgrid &&
       roundPatches &&
       spacing(recoGrid) == spacing(grids[k])

      issubgrid = true
      grids[k] = deriveSubgrid(recoGrid, grids[k])

    end

    # if the patch is a true subgrid we don't need to apply interpolation and can load the
    # matrix as is.
    if issubgrid
      # we first check if the matrix is already in memory
      u = -1
      for l=1:length(SOrigIdx)
        if SOrigIdx[l] == idx && SIsPlain[l]
          u = l
          break
        end
      end
      if u > 0 # its already in memory
        patchToSMIdx[k] = u
      else     # not yet in memory  -> load it
        S_, grid = getSF(SF,freq,nothing,"Kaczmarz", bgCorrection=bgCorrection, tfCorrection=tfCorrection)
        push!(S,S_)
        push!(SOrigIdx,idx)
        push!(SIsPlain,true) # mark this as a plain system matrix (without interpolation)
        patchToSMIdx[k] = length(S)
      end
    else
      # in this case the patch grid does not fit onto the reco grid. Lets derive a subgrid
      # that is very similar to grids[k]
      newGrid = deriveSubgrid(recoGrid, grids[k])

      # load the matrix on the new subgrid
      S_, grid = getSF(SF,freq,nothing,"Kaczmarz", bgCorrection=bgCorrection,
                   gridsize=shape(newGrid),
                   fov=fieldOfView(newGrid),
                   center=fieldOfViewCenter(newGrid).-fieldOfViewCenter(grids[k]),
		   tfCorrection=tfCorrection)
                   # @TODO: I don't know the sign of aboves statement

      grids[k] = newGrid # we need to change the stored Grid since we now have a true subgrid
      push!(S,S_)
      push!(SOrigIdx,idx)
      push!(SIsPlain,false)
      patchToSMIdx[k] = length(S)
    end
  end
  @debug "Use $(length(S)) patches"

  @debug "Calculate LUT"
  # now that we have all grids we can calculate the indices within the recoGrid
  xcc, xss = calculateLUT(grids, recoGrid)

  sign = ones(Int, M, numPatches)

  return MultiPatchOperator{eltype(first(S)), reduce(promote_type, typeof.(S)), typeof(recoGrid), eltype(patchToSMIdx), typeof(patchToSMIdx), typeof(sign)}(S, recoGrid, length(recoGrid), M*numPatches,
             RowToPatch, xcc, xss, sign, numPatches, patchToSMIdx)
end


function calculateLUT(grids, recoGrid)
  xss = Vector{Int}[]
  xcc = Vector{Int}[]
  for k=1:length(grids)
    N = length(grids[k])
    push!(xss, collect(1:N))
    xc = zeros(Int64,N)
    for n=1:N
      xc[n] = posToLinIdx(recoGrid,grids[k][n])
    end
    push!(xcc, xc)
  end
  return xcc, xss
end

function size(FFOp::MultiPatchOperator,i::Int)
  if i==2
    return FFOp.N
  elseif i==1
    return FFOp.M
  else
    error("bounds error")
  end
end

size(FFOp::MultiPatchOperator) = (FFOp.M,FFOp.N)
size(FFTOp::DenseMultiPatchOperator) = (FFTOp.M,FFTOp.N)

length(FFOp::MultiPatchOperator) = size(FFOp,1)*size(FFOp,2)
function getindex(op::MultiPatchOperator, i::I, j::I) where I <: Integer
  p = op.RowToPatch[i]
  xs = op.xss[p]
  xc = op.xcc[p]
  row = mod1(i,div(op.M,op.nPatches))
  A = op.S[op.patchToSMIdx[p]]
  sign = op.sign[row,op.patchToSMIdx[p]]
  index = findfirst(isequal(j), xc)
  if !isnothing(index)
    return sign*A[row,xs[index]]
  else
    return zero(eltype(op))
  end
end

function LinearAlgebra.mul!(b::AbstractVector{T}, op::MultiPatchOperator{T}, x::AbstractVector{T}) where T
  for i in 1:size(op, 1)
    b[i] = dot_with_matrix_row(op, x, i)
  end
  return b
end

function LinearAlgebra.mul!(res::AbstractVector{T}, adj::Adjoint{T, OP}, t::AbstractVector{T}) where {T, V <: AbstractArray, OP <: MultiPatchOperator{T, V}}
  op = adj.parent
  res .= zero(T)

  for i in 1:size(op, 1)
    val = t[i]

    p = op.RowToPatch[i]
    xs = op.xss[p]
    xc = op.xcc[p]
    row = mod1(i,div(op.M,op.nPatches))
    A = op.S[op.patchToSMIdx[p]]
    sign = op.sign[row,op.patchToSMIdx[p]]
  
    for j in 1:length(xs)
      res[xc[j]] += adjoint(sign*A[row,xs[j]]) * val
    end
  end
  return res
end

### The following is intended to use the standard kaczmarz method ###
function RegularizedLeastSquares.normalize(norm::SystemMatrixBasedNormalization, op::MultiPatchOperator, b)
  if length(op.S) == 1
    trace = RegularizedLeastSquares.normalize(norm, op.S[1], b)
    trace *= op.nPatches #*prod(Op.PixelSizeSF)/prod(Op.PixelSizeC)
  else
    trace = sum([RegularizedLeastSquares.normalize(norm, S, b)*size(S, 2) for S in op.S])
    #trace *= prod(Op.PixelSizeSF)/prod(Op.PixelSizeC)
    trace/=size(op, 2)
  end
  return trace
end

function RegularizedLeastSquares.normalize(norm::SystemMatrixBasedNormalization, op::DenseMultiPatchOperator, b)
  trace = sum([RegularizedLeastSquares.normalize(norm, view(op.S, :, :, i), b)*size(op.S, 2) for i in axes(op.S, 3)])
  trace/=size(op, 2)
  return trace
end

function calculateTraceOfNormalMatrix(Op::MultiPatchOperator, weights)
  if length(Op.S) == 1
    trace = calculateTraceOfNormalMatrix(Op.S[1],weights)
    trace *= Op.nPatches #*prod(Op.PixelSizeSF)/prod(Op.PixelSizeC)
  else
    trace = sum([calculateTraceOfNormalMatrix(S,weights) for S in Op.S])
    #trace *= prod(Op.PixelSizeSF)/prod(Op.PixelSizeC)
  end
  return trace
end

setlambda(::MultiPatchOperator, ::Any) = nothing

function dot_with_matrix_row(Op::MultiPatchOperator, x::AbstractArray{T}, k::Integer) where T
  p = Op.RowToPatch[k]
  xs = Op.xss[p]
  xc = Op.xcc[p]

  j = mod1(k,div(Op.M,Op.nPatches))
  A = Op.S[Op.patchToSMIdx[p]]
  sign = Op.sign[j,Op.patchToSMIdx[p]]

  return dot_with_matrix_row_(A,x,xs,xc,j,sign)
end

function dot_with_matrix_row_(A::AbstractArray{T},x,xs,xc,j,sign) where T
  tmp = zero(T)
  @simd  for i = 1:length(xs)
     @inbounds tmp += sign*A[j,xs[i]]*x[xc[i]]
  end
  tmp
end

function kaczmarz_update!(Op::MultiPatchOperator, x::AbstractArray, k::Integer, beta)
  p = Op.RowToPatch[k]
  xs = Op.xss[p]
  xc = Op.xcc[p]

  j = mod1(k,div(Op.M,Op.nPatches))
  A = Op.S[Op.patchToSMIdx[p]]
  sign = Op.sign[j,Op.patchToSMIdx[p]]

  kaczmarz_update_!(A,x,beta,xs,xc,j,sign)
end

function kaczmarz_update_!(A,x,beta,xs,xc,j,sign)
  @simd for i = 1:length(xs)
    @inbounds x[xc[i]] += beta* conj(sign*A[j,xs[i]])
  end
end

function RegularizedLeastSquares.rownorm²(op::MultiPatchOperator, row::Int64)
  p = op.RowToPatch[row]
  xs = op.xss[p]

  j = mod1(row,div(op.M,op.nPatches))
  A = op.S[op.patchToSMIdx[p]]
  sign = op.sign[j,op.patchToSMIdx[p]]

  return mapreduce(x -> abs2(sign*A[j,x]), +, xs)
end
