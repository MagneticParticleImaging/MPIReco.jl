import Base: size
import LinearSolver: initkaczmarz, dot_with_matrix_row, kaczmarz_update!

export reconstructionSeparate, reconstructionFFJoint, FFOperator

# this wrapper should maybe be removed
function reconstruction{T<:MPIFile}(bSF::MPIFile, bMeas::Vector{T}, freq; kargs...)
  reconstructionFFJoint(bSF, bMeas, freq; kargs...)
end

# Necessary for Multi-System-Matrix FF reconstruction
voxelSize(bSF::Vector{T}) where T<:MPIFile = voxelSize(bSF[1])
sfGradient(bSF::Vector{T},dim) where T<:MPIFile = gradient(bSF[1],dim)
generateHeaderDict(bSF::Vector{T},bMeas::MPIFile) where T<:MPIFile =
   generateHeaderDict(bSF[1],bMeas)

function reconstructionFFJoint(bSF, bMeas::MPIFile, freq;
            frames=nothing, bEmpty=nothing, OverscanSF=[0,0,0],OffsetFF=[0.002,0.002,0.002],
            nAverages=1,denoiseWeight=0, loadas32bit=false, alpha=[0,0,0],
            spectralLeakageCorrection=true, kargs...)

  consistenceCheck(bSF, bMeas)

  bgcorrection = (bEmpty != nothing)

  FFPos = ffPos(bMeas,alpha=alpha)

  periodsSortedbyFFPos = unflattenOffsetFieldShift(FFPos)

  FFPos = FFPos[:,periodsSortedbyFFPos[:,1]]

  FFOp = FFOperator(bSF,bMeas,freq,bgcorrection,OverscanSF=OverscanSF,
                    OffsetFF=OffsetFF,denoiseWeight=denoiseWeight,FFPos=FFPos,indFFPos=periodsSortedbyFFPos[:,1])

  L = numScans(bMeas)
  (frames==nothing) && (frames=collect(1:L))
  nFrames=length(frames)

  uTotal = getMeasurementsFD(bMeas,frequencies=freq, frames=frames, numAverages=nAverages,
                             spectralLeakageCorrection=spectralLeakageCorrection)

  uTotal=uTotal[:,periodsSortedbyFFPos,:]
  uTotal=mean(uTotal,3)

  # Here we call a regular reconstruction function
  c = reconstruction(FFOp,uTotal,(FFOp.PixelSizeC...,); kargs...)

  pixspacing = voxelSize(bSF) ./ sfGradient(bMeas,3) .* sfGradient(bSF,3)
  offset = mean(FFOp.CSize_mm, 2)  .- 0.5.*pixspacing.*FFOp.PixelSizeC .+ 0.5.*pixspacing

  im = AxisArray(c, (:x,:y,:z,:time), tuple(pixspacing...,dfcycle(bMeas)),
                                      tuple(offset...,0.0))

  imMeta = ImageMeta(im,generateHeaderDict(bSF,bMeas))
  return imMeta
end


# FFOperator is a type that acts as the MPI system matrix but exploits
# its sparse structure.
# Its very important to keep this type typestable
type FFOperator{V<:AbstractMatrix}
  S::Vector{V}
  PixelSizeC::Vector{Int}
  PixelSizeSF::Vector{Int}
  SFIndex::Matrix{Int}
  CIndex::Matrix{Int}
  Length::Matrix{Int}
  Listc::Vector{Int}
  Lists::Vector{Int}
  CSize_mm::Matrix{Float64}
  freq::Vector{Int}
  N::Int
  M::Int
  ProdL::Matrix{Int}
  RowToPatch::Vector{Int}
  xcc::Vector{Vector{Int}}
  xss::Vector{Vector{Int}}
  nPatches::Int
  patchToSMIdx::Vector{Int}
end

function FFOperator(SF::MPIFile, bMeas,freq::Vector{Int64},bgcorrection::Bool;
                    OverscanSF=nothing, OffsetFF=nothing, denoiseWeight=0,
                    FFPos=ffPos(bMeas),indFFPos=nothing, kargs...)

   indFFPos == nothing ? indFFPos=analyseFFPos(FFPos) : nothing

   S = getSF(SF,freq,nothing,"kaczmarz", bgcorrection=bgcorrection)

   shape = getshape(gridSize(SF))

   if denoiseWeight > 0
     denoiseSF!(S, shape, weight=denoiseWeight)
   end

   CIndex=zeros(Int64,length(indFFPos),3)
   SFIndex=zeros(Int64,length(indFFPos),3)

   Length=zeros(Int64,length(indFFPos),3)
   Listc=Int64[]
   Lists=Int64[]

   PixelSizeC=zeros(Int64,3)

   voxelS=voxelSize(SF)
   PixelSizeSF=round.(Int64,fov(SF)./voxelS)

   OverscanSF==nothing ? SFSize_mm=SFSize(SF) : SFSize_mm=SFSize(SF,OverscanSF=OverscanSF)
   OffsetFF==nothing ? CSize_mm=MultiPatchImageSize(dfFov(bMeas),FFPos) : CSize_mm=MultiPatchImageSize(dfFov(bMeas),FFPos,OffsetFF=OffsetFF)
   rr=8

   fovCenter = FFPos#(bMeas,alpha=alpha)
   for i = 1:length(indFFPos)#numPatches(bMeas)
     areaSF=zeros(Float64,3,2)
     areaObj=zeros(Float64,3,2)
     for k=1:3

        #PixelSizeC[k]=ceil(Int,CSize_mm[k,2]-CSize_mm[k,1])/voxelS[k])
        PixelSizeC[k]=round.(Int, (CSize_mm[k,2]-CSize_mm[k,1])/voxelS[k])

        PixSF=round.(linspace(-fov(SF)[k]/2+voxelS[k]/2,fov(SF)[k]/2-voxelS[k]/2,round.(Int32,fov(SF)[k]/voxelS[k])),rr)

        PixC=linspace(CSize_mm[k,1],CSize_mm[k,2],PixelSizeC[k]+1)

        areaSF[k,:]=round.([maximum([-SFSize_mm[k]/2,CSize_mm[k,1]-fovCenter[k,i]]);minimum([SFSize_mm[k]/2,CSize_mm[k,2]-fovCenter[k,i]])],rr)
        areaObj[k,:]=round.(areaSF[k,:]+fovCenter[k,i],rr)
        help=collect(areaObj[k,1]:voxelS[k]:areaObj[k,2])

        tmp1=Int32[]
        tmp2=Int32[]
        tmp3=Int32[]
        for j=1:maximum([length(PixSF),length(PixC)])
            (j<=length(PixSF) && PixSF[j]>=areaSF[k,1] && PixSF[j]<=areaSF[k,2]) ? push!(tmp1,j):nothing
            (j<=length(PixC) && PixC[j]>=areaObj[k,1] && PixC[j]<=areaObj[k,2]) ? push!(tmp2,j):nothing
        end

        for p=1:PixelSizeC[k]
            for l=1:length(help)
               abs(PixC[p]-help[l])<voxelS[k]/2 ? push!(tmp3,p) : nothing
            end
        end
       tmp2=tmp3
       #k==1 && tmp1[1]!=1 ? tmp1=tmp1-1 : nothing

       #k==3 && tmp1[end]!=fov(SF)[3]/voxelS[3] ? tmp1=tmp1+1 : nothing
        tmp1!=[] ? SFIndex[i,k]=tmp1[1] : SFIndex[i,k]=0

        tmp2!=[] ? CIndex[i,k]=tmp2[1] : CIndex[i,k]=0

        Length[i,k]=minimum([length(tmp1),length(tmp2)])
    end
    xc=floor(Int, CIndex[i,1]+PixelSizeC[1]*(CIndex[i,2]-1)+PixelSizeC[1]*PixelSizeC[2]*(CIndex[i,3]-1))

    xs=floor(Int, SFIndex[i,1]+PixelSizeSF[1]*(SFIndex[i,2]-1)+PixelSizeSF[1]*PixelSizeSF[2]*(SFIndex[i,3]-1))

    kronic=collect(kron(ones(Int64,Length[i,2]*Length[i,3]),collect(0:Length[i,1]-1)).+PixelSizeC[1]*kron(ones(Int64,Length[i,3]),collect(0:Length[i,2]-1),ones(Int64,Length[i,1])).+PixelSizeC[1]*PixelSizeC[2]*kron(collect(0:Length[i,3]-1),ones(Int64,Length[i,2]*Length[i,1])))
    kronis=collect(kron(ones(Int64,Length[i,2]*Length[i,3]),collect(0:Length[i,1]-1)).+PixelSizeSF[1]*kron(ones(Int64,Length[i,3]),collect(0:Length[i,2]-1),ones(Int64,Length[i,1])).+PixelSizeSF[1]*PixelSizeSF[2]*kron(collect(0:Length[i,3]-1),ones(Int64,Length[i,2]*Length[i,1])))
    Listc=[Listc;kronic+xc]

    Lists=[Lists;kronis+xs]

   end

  N=prod(PixelSizeSF)
  M=div(length(S),N)

  numPatches = size(SFIndex)[1]
  RowToPatch = kron(collect(1:numPatches), ones(Int,M))

  ProdL = [0;prod(Length,2)]

  xss = Vector{Int}[]
  xcc = Vector{Int}[]
  for k=1:size(SFIndex)[1]
    push!(xcc, Listc[1+sum(ProdL[1:k]):sum(ProdL[1:k+1])])#Op.Lists
    push!(xss, Lists[1+sum(ProdL[1:k]):sum(ProdL[1:k+1])])#Op.Listc
  end


  patchToSMIdx = ones(Int, numPatches)

  FFOperator([S],PixelSizeC,PixelSizeSF, SFIndex, CIndex, Length, Listc, Lists, CSize_mm,freq,
      prod(PixelSizeC), M*(size(SFIndex)[1]), ProdL, RowToPatch, xcc, xss, numPatches, patchToSMIdx)
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



function MultiPatchImageSize(DFFOV,FFPos; OffsetFF=[0.0,0.0,0.0])
   nP=size(FFPos,2)
   ffPositions = zeros(nP*2,3)
   ResultingSize= zeros(3,2)

   for i=1:nP
      if DFFOV[:,1] != DFFOV[:,i]
        warn("DF FoV is not equal")
      end
      ffPositions[i,:] = FFPos[:,i] + DFFOV[:,1]./2
      ffPositions[i+nP,:] = FFPos[:,i] - DFFOV[:,1]./2
   end
   for k=1:3
      ResultingSize[k,:]=round.([(minimum(ffPositions[:,k]))-OffsetFF[k],((maximum(ffPositions[:,k]))+OffsetFF[k])],8)
   end
    return ResultingSize
end

# Offset to the dfFov, Nothing means the whole SFSize
function SFSize(SF::BrukerFile; OverscanSF = [0.0,0.0,0.0])
    if OverscanSF == [0,0,0]
       OverscanSF=(fov(SF)-dfFov(SF))./2
    else
       length(OverscanSF)<3 ? OverscanSF=[OverscanSF,zeros(3-length(OverscanSF))]: OverscanSF=OverscanSF[1:3]
    end

    if sum(2*OverscanSF+dfFov(SF).>fov(SF))!=0
        println("FOV=",fov(SF),", ","DF=",dfFov(SF)," and ","DF+Offset=",2*OverscanSF+dfFov(SF))
        error("To Big SF Offset!")
    end

    return round.(dfFov(SF)+2*OverscanSF,8)
end

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
        s² = rownorm²(Op.S[l],i)*weights[i]^2
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
