import Base.length, Base.size

export getSF, SVD, tikhonovLU, setlambda

include("SystemMatrixRecovery.jl")
include("SystemMatrixCenter.jl")

function converttoreal(S::AbstractArray{Complex{T}},f) where T
  N = prod(calibSize(f))
  M = div(length(S),N)
  S = reshape(S,N,M)
  S = reshape(reinterpret(T,vec(S)),(2*N,M))
  for l=1:M
    tmp = S[:,l]
    S[1:N,l] = tmp[1:2:end]
    S[N+1:end,l] = tmp[2:2:end]
  end
  return S
end

setlambda(S::AbstractMatrix, λ) = nothing

function getSF(bSF, frequencies, sparseTrafo, solver; kargs...)
  SF, grid = getSF(bSF, frequencies, sparseTrafo; kargs...)
  if solver == "kaczmarz"
    return transpose(SF), grid
  elseif solver == "pseudoinverse"
    return SVD(svd(transpose(SF))...), grid
  elseif solver == "cgnr" || solver == "lsqr" || solver == "fusedlasso"
    return copy(transpose(SF)), grid
  elseif solver == "direct"
    return RegularizedLeastSquares.tikhonovLU(copy(transpose(SF))), grid
  else
    return SF, grid
  end
end

getSF(bSF::Union{T,Vector{T}}, frequencies, sparseTrafo::Nothing; kargs...) where {T<:MPIFile} =
   getSF(bSF, frequencies; kargs...)

# TK: The following is a hack since colored and sparse trafo are currently not supported
getSF(bSF::Vector{T}, frequencies, sparseTrafo::AbstractString; kargs...) where {T<:MPIFile} =
  getSF(bSF[1], frequencies, sparseTrafo; kargs...)

getSF(bSFs::MultiContrastFile, frequencies; kargs...) =
   getSF(bSFs.files, frequencies; kargs...)

function getSF(bSFs::Vector{T}, frequencies; kargs...) where {T<:MPIFile}
  maxFov = [0.0,0.0,0.0]
  maxSize = [0,0,0]
  for l=1:length(bSFs)
    maxFov = max.(maxFov, calibFov(bSFs[l]))
    maxSize = max.(maxSize, collect(calibSize(bSFs[l])) )
  end
  data = [getSF(bSF, frequencies; gridsize=maxSize, fov=maxFov, kargs...) for bSF in bSFs]
  return cat([d[1] for d in data]..., dims=1), data[1][2] # grid of the first SF
end

getSF(f::MPIFile; recChannels=1:numReceivers(f), kargs...) = getSF(f, filterFrequencies(f, recChannels=recChannels); kargs...)

#=
function repairDeadPixels(S, shape, deadPixels)
  shapeT = tuple(shape...)
  for k=1:size(S,2)
    i2s = CartesianIndices(shapeT)
    for dp in deadPixels
      ix,iy,iz=i2s[dp].I
      if 1<ix<shape[1] && 1<iy<shape[2] && 1<iz<shape[3]
         lval = S[ LinearIndices(shapeT)[ix,iy+1,iz]  ,k]
         rval = S[ LinearIndices(shapeT)[ix,iy-1,iz]  ,k]
         S[dp,k] = 0.5*(lval+rval)
         fval = S[ LinearIndices(shapeT)[ix+1,iy,iz]  ,k]
         bval = S[ LinearIndices(shapeT)[ix-1,iy,iz]  ,k]
         S[dp,k] = 0.5*(fval+bval)
         uval = S[ LinearIndices(shapeT)[ix,iy,iz+1]  ,k]
         dval = S[ LinearIndices(shapeT)[ix,iy,iz-1]  ,k]
         S[dp,k] = 0.5*(uval+dval)
      end
    end
  end
end
=#
function getSF(bSF::MPIFile,saveto::String;recChannels=1:numReceivers(bSF),kargs...)
  S,grid = getSF(bSF,filterFrequencies(bSF, recChannels=recChannels);bgCorrection=false,returnasmatrix=true,kargs...)
  if !isempty(saveto)
    @info "Saving transformed SystemMatrix..."
    params=MPIFiles.loadDataset(bSF)
    bG = params[:measData][acqNumFGFrames(bSF)+1:end,:,:,:]
    params[:measData]=vcat(reshape(S,size(S,1),Int(size(S,2)/3),3,1),bG)
    params[:acqNumFrames]=size(vcat(reshape(S,size(S,1),Int(size(S,2)/3),3,1),bG),1)
    params[:calibFov]=grid.fov
    params[:calibSize]=grid.shape
    params[:calibFovCenter]=grid.center
    # Wann ist diese Zuordnung nicht einfach fortlaufend..dieser Fall ist so nicht abgedeckt
    params[:measIsBGFrame]=(i-> i in collect(size(S,1)+1:size(S,1)+size(bG,1))).(collect(1:size(S,1)+size(bG,1)))
    saveasMDF(saveto, params)
  end
  return S,grid
end

function getSF(bSF::MPIFile, frequencies; returnasmatrix = true, procno::Integer=1,
               bgcorrection=false, bgCorrection=bgcorrection, loadasreal=false,
	             gridsize=collect(calibSize(bSF)), numPeriodAverages=1,numPeriodGrouping=1,
	             fov=calibFov(bSF), center=[0.0,0.0,0.0],deadPixels=nothing, kargs...)

  nFreq = rxNumFrequencies(bSF)

  numPeriods = div(acqNumPeriodsPerFrame(bSF),numPeriodGrouping*numPeriodAverages)

  S = getSystemMatrix(bSF, frequencies, bgCorrection=bgCorrection,
                      numPeriodAverages=numPeriodAverages,
                      numPeriodGrouping=numPeriodGrouping; kargs...)

  if deadPixels != nothing
    #repairDeadPixels(S,gridsize,deadPixels)
    @info "Repairing deadPixels..."
    S = repairSM(S,RegularGridPositions(calibSize(bSF),calibFov(bSF),[0.0,0.0,0.0]),deadPixels)
  end

  if collect(gridsize) != collect(calibSize(bSF)) ||
    center != [0.0,0.0,0.0] ||
    fov != calibFov(bSF)

    origin = RegularGridPositions(calibSize(bSF),calibFov(bSF),[0.0,0.0,0.0])
    target = RegularGridPositions(gridsize,fov,center)

    if any(fov .> calibFov(bSF))
      #round.(Int,(fov .- origin.fov).*(origin.shape./(2 .* origin.fov)))
      if gridsize == collect(calibSize(bSF))
        gridsize_new = round.(Int, fov .* origin.shape ./ (2 * origin.fov)) * 2
        @info "You selected a customized (bigger) FOV, without selecting a bigger grid. Thus, an Extrapolation to
the new FOV is followed by a Interpolation to the old grid-size, leading to a change in gridpoint-size. If you want
to roughly keep the original gridpoint-size, define the key-word gridsize = $gridsize_new,
alongside to your FOV-selection."
      end
      S,origin = extrapolateSM(S,origin,fov)
      # alternativ ginge direkt: extrapolateSM(SM, grid, fov)
    end

    @info "Perform SF Interpolation..."

    SInterp = zeros(eltype(S),prod(gridsize),length(frequencies)*numPeriods)
    for k=1:length(frequencies)*numPeriods
      A = MPIFiles.interpolate(reshape(S[:,k],origin.shape...), origin, target)
      SInterp[:,k] = vec(A)
    end
    S = SInterp
    grid = target
  else
    grid = RegularGridPositions(calibSize(bSF),calibFov(bSF),[0.0,0.0,0.0])
  end

  if loadasreal
    S = converttoreal(S,bSF)
    resSize = [gridsize..., 2*length(frequencies)*numPeriods]
  else
    resSize = [gridsize...,length(frequencies)*numPeriods]
  end

  if returnasmatrix
    return reshape(S,prod(resSize[1:end-1]),resSize[end]), grid
  else
    return squeeze(reshape(S,resSize...)), grid
  end
end

"""
This function reads the component belonging to the given mixing-factors und channel from the given MPIfile.\n
In the case of a 2D systemfunction, the third dimension is added.
"""
function getSF(bSF::MPIFile, mx::Int, my::Int, mz::Int, recChan::Int; bgcorrection=true)
  maxFreq  = rxNumFrequencies(bSF)
  freq = mixFactorToFreq(bSF, mx, my, mz)
  freq = clamp(freq,0,maxFreq-1)
  k = freq + (recChan-1)*maxFreq
  @debug "Frequency = $k"
  SF, grid = getSF(bSF,[k+1],returnasmatrix = true, bgcorrection=bgcorrection)
  N = calibSize(bSF)
  SF = reshape(SF,N[1],N[2],N[3])

  return SF
end

include("Sparse.jl")
