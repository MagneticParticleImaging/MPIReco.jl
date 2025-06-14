import Base.length, Base.size

export getSF, SVD, tikhonovLU, setlambda

include("SystemMatrixRecovery.jl")
include("SystemMatrixCenter.jl")
include("SystemMatrixWrapper.jl")
include("SMExtrapolation.jl")

export AbstractSystemMatrixParameter
abstract type AbstractSystemMatrixParameter <: AbstractMPIRecoParameters end

export AbstractSystemMatrixGriddingParameter
abstract type AbstractSystemMatrixGriddingParameter <: AbstractSystemMatrixParameter end

export NoGridding
struct NoGridding <: AbstractSystemMatrixGriddingParameter end

export SystemMatrixGriddingParameter
Base.@kwdef struct SystemMatrixGriddingParameter <: AbstractSystemMatrixGriddingParameter
  gridsize::Vector{Int64} = [1, 1, 1]
  fov::Vector{Float64} = [0.0, 0.0, 0.0]
  center::Vector{Float64} = [0.0,0.0,0.0]
  deadPixels::Union{Nothing, Vector{Int64}} = nothing
end

SystemMatrixGriddingParameter(file::MPIFile) = SystemMatrixGriddingParameter(
  gridsize = calibSize(file),
  fov = calibFov(file),
  center = calibFovCenter(file)
)

export defaultParameterTfCorrection, defaultParameterGridSize, defaultParameterCalibCenter, defaultParameterCalibFov
defaultParameterTfCorrection(new::MPIFile) = rxHasTransferFunction(new)
defaultParameterTfCorrection(new::Missing) = missing
defaultParameterGridSize(new::MPIFile) = gridSizeCommon(new)
defaultParameterGridSize(new::Missing) = missing
defaultParameterCalibCenter(new::MPIFile) = calibFovCenter(new)
defaultParameterCalibCenter(new::Missing) = missing
defaultParameterCalibFov(new::MPIFile) = calibFov(new)
defaultParameterCalibFov(new::Missing) = missing

# Maybe implement custom defaults with optional given sf -> remove @kwdef
#function SystemMatrixGriddingParameter(;sf::MPIFile, gridsize = nothing, fov = nothing, center = [0.0, 0.0, 0.0], deadPixels = Int64[])
#  if isnothing(gridsize)
#    gridsize = gridSizeCommon(sf)
#  end
#  if isnothing(fov)
#    fov = calibFov(sf)
#  end
#  return SNRThresholdFrequencyFilterParameter(gridsize, fov, center, deadPixels)
#end
export AbstractSystemMatrixLoadingParameter
abstract type AbstractSystemMatrixLoadingParameter <: AbstractSystemMatrixParameter end

export DenseSystemMatixLoadingParameter
Base.@kwdef struct DenseSystemMatixLoadingParameter{F<:AbstractFrequencyFilterParameter, G<:AbstractSystemMatrixGriddingParameter} <: AbstractSystemMatrixLoadingParameter
  freqFilter::F
  gridding::G
  bgCorrection::Bool = false
  tfCorrection::Union{Bool, Nothing} = nothing
  loadasreal::Bool = false
end
function process(t::Type{<:AbstractMPIRecoAlgorithm}, params::DenseSystemMatixLoadingParameter, sf::MPIFile)
  # Construct freqFilter
  frequencies = process(t, params.freqFilter, sf)
  return frequencies, process(t, params, sf, frequencies)...
end
function process(t::Type{<:AbstractMPIRecoAlgorithm}, params::DenseSystemMatixLoadingParameter, sf::MPIFile, frequencies::Vector{CartesianIndex{2}})
  S, grid = getSF(sf, frequencies, nothing; toKwargs(params, default = Dict{Symbol, Any}(:tfCorrection => rxHasTransferFunction(sf)))...)
  @info "Loading SM"
  return S, grid
end

export SparseSystemMatrixLoadingParameter
Base.@kwdef struct SparseSystemMatrixLoadingParameter{F<:AbstractFrequencyFilterParameter} <: AbstractSystemMatrixLoadingParameter
  freqFilter::F
  sparseTrafo::Union{Nothing, String} = nothing # TODO concrete options here?
  thresh::Float64 = 0.0
  redFactor::Float64 = 0.1
  bgCorrection::Bool = false
  tfCorrection::Union{Bool, Nothing} = nothing
  loadasreal::Bool=false
  useDFFoV::Bool = false
end
function process(t::Type{<:AbstractMPIRecoAlgorithm}, params::SparseSystemMatrixLoadingParameter, sf::MPIFile)
  # Construct freqFilter
  frequencies = process(t, params.freqFilter, sf)
  return frequencies, process(t, params, sf, frequencies)...
end
function process(t::Type{<:AbstractMPIRecoAlgorithm}, params::SparseSystemMatrixLoadingParameter, sf::MPIFile, frequencies::Vector{CartesianIndex{2}})
  S, grid = getSF(sf, frequencies, params.sparseTrafo; toKwargs(params, default = Dict{Symbol, Any}(:tfCorrection => rxHasTransferFunction(sf)))...)
  return S, grid
end
function process(t::Type{<:AbstractMPIRecoAlgorithm}, params::SparseSystemMatrixLoadingParameter, elType::Type{<:Number}, arrayType, shape::NTuple{N, Int64}) where N
  return createLinearOperator(params.sparseTrafo, elType; shape, S = typeof(arrayType{elType}(undef, 0)))
end

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

function getSF(bSF, frequencies, sparseTrafo, solver::AbstractString; kargs...)
  if solver == "Kaczmarz"
    return getSF(bSF, frequencies, sparseTrafo, Kaczmarz; kargs...)
  elseif solver == "pseudoinverse"
    return getSF(bSF, frequencies, sparseTrafo, PseudoInverse; kargs...)
  elseif solver == "cgnr" || solver == "fusedlasso"
    return getSF(bSF, frequencies, sparseTrafo, CGNR; kargs...)
  elseif solver == "direct"
    return getSF(bSF, frequencies, sparseTrafo, DirectSolver; kargs...)
  else
    @warn "Unsupported solver $(solver) used in getSF"
    return getSF(bSF, frequencies, sparseTrafo; kargs...)
  end
end
getSF(bSF, frequencies, sparseTrafo, solver::AbstractLinearSolver; kargs...) = getSF(bSF, frequencies, sparseTrafo, typeof(solver); kargs...)
function getSF(bSF, frequencies, sparseTrafo, solver::Type{<:AbstractLinearSolver}; arrayType = Array, kargs...)
  SF, grid = getSF(bSF, frequencies, sparseTrafo; kargs...)
  SF, grid = prepareSF(solver, SF, grid)
  SF = adaptSF(arrayType, SF)
  return SF, grid
end

# In this instance we want to dispatch before the cache and call individual steps processing steps which in turn can use the cache. This is only effective if the cache size is >= number of sub-processes
# Alternative solutions require us to nest our parameters with caches in between
function AbstractImageReconstruction.process(type::Type{<:AbstractMPIRecoAlgorithm}, params::Union{L, ProcessResultCache{L}}, sf::MPIFile, solverT, arrayType = Array) where L <: AbstractSystemMatrixLoadingParameter
  # Each process step can access the cache
  freqs, sf, grid = process(type, params, sf)
  sf, grid = process(type, params, sf, solverT, grid)
  sf = process(type, params, sf, arrayType)
  return freqs, sf, grid
end
function AbstractImageReconstruction.process(type::Type{<:AbstractMPIRecoAlgorithm}, params::L, sf::AbstractArray, solverT::Type{<:AbstractLinearSolver}, grid) where L <: AbstractSystemMatrixLoadingParameter
  @info "Preparing SF"
  return prepareSF(solverT, sf, grid)
end
function AbstractImageReconstruction.process(type::Type{<:AbstractMPIRecoAlgorithm}, params::L, sf::AbstractArray, arrayType::Type{<:AbstractArray}) where L <: AbstractSystemMatrixLoadingParameter
  @info "Adapting SF"
  return adaptSF(arrayType, sf)
end

# Assumption SF is a (wrapped) CPU-array
# adapt(Array, Sparse-CPU) results in dense array, so we only want to adapt if necessary
adaptSF(arrayType, SF) = adapt(arrayType, SF)
adaptSF(arrayType::Type{<:Array}, SF) = SF

prepareSF(solver::Type{Kaczmarz}, SF, grid) = transpose(SF), grid
prepareSF(solver::Type{PseudoInverse}, SF, grid) = SVD(svd(transpose(SF))...), grid
prepareSF(solver::Type{DirectSolver}, SF, grid) = RegularizedLeastSquares.tikhonovLU(copy(transpose(SF))), grid
prepareSF(solver::Type{Kaczmarz}, SF::AbstractSparseArray, grid) = transpose(SF), grid
prepareSF(solver::Type{PseudoInverse}, SF::AbstractSparseArray, grid) = SVD(svd(transpose(SF))...), grid
prepareSF(solver::Type{DirectSolver}, SF::AbstractSparseArray, grid) = RegularizedLeastSquares.tikhonovLU(copy(transpose(SF))), grid

prepareSF(solver, SF , grid) = copy(transpose(SF)), grid
prepareSF(solver, SF::AbstractSparseArray, grid) = sparse(transpose(SF)), grid


prepareNormalSF(solver::AbstractLinearSolver, SF) = prepareNormalSF(typeof(solver), SF)
prepareNormalSF(solver::Type{<:RegularizedLeastSquares.AbstractLinearSolver}, SF) = nothing
prepareNormalSF(solver::Union{Type{FISTA}, Type{OptISTA}, Type{POGM}, Type{CGNR}, Type{ADMM}, Type{SplitBregman}}, SF) = LinearOperatorCollection.normalOperator(SF)

function getSF(bSF::Union{T,Vector{T}}, frequencies, sparseTrafo::Nothing; kargs...) where {T<:MPIFile}
  return getSF(bSF, frequencies; kargs...)
end

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
	             fov=calibFov(bSF), center=calibFovCenter(bSF),deadPixels=nothing, kargs...)

  nFreq = rxNumFrequencies(bSF)

  numPeriods = div(acqNumPeriodsPerFrame(bSF),numPeriodGrouping*numPeriodAverages)

  S = getSystemMatrix(bSF, frequencies, bgCorrection=bgCorrection,
                      numPeriodAverages=numPeriodAverages,
                      numPeriodGrouping=numPeriodGrouping; kargs...)

  if deadPixels != nothing
    #repairDeadPixels(S,gridsize,deadPixels)
    @info "Repairing deadPixels..."
    S = repairSM(S,RegularGridPositions(calibSize(bSF),calibFov(bSF),calibFovCenter(bSF)),deadPixels)
  end

  if collect(gridsize) != collect(calibSize(bSF)) ||
    center != calibFovCenter(bSF) ||
    fov != calibFov(bSF)

    origin = RegularGridPositions(calibSize(bSF),calibFov(bSF),calibFovCenter(bSF))
    target = RegularGridPositions(gridsize,fov,center)

    if any(fov .> calibFov(bSF))
      #round.(Int,(fov .- origin.fov).*(origin.shape./(2 .* origin.fov)))
      if gridsize == collect(calibSize(bSF))
        gridsize_new = round.(Int, fov .* origin.shape ./ (2 * origin.fov),RoundNearestTiesUp) * 2
        @info "You selected a customized (bigger) FOV, without selecting a bigger grid. Thus, an Extrapolation to
the new FOV is followed by an Interpolation to the old grid-size, leading to a change in gridpoint-size. If you want
to roughly keep the original gridpoint-size, define the key-word gridsize = $gridsize_new,
alongside to your FOV-selection."
      end
      S,origin = extrapolateSM(S,origin,fov)
      # alternativ ginge direkt: extrapolateSM(SM, grid, fov)
    end

    @debug "Perform SF Interpolation..."

    SInterp = zeros(eltype(S),prod(gridsize),length(frequencies)*numPeriods)
    for k=1:length(frequencies)*numPeriods
      A = MPIFiles.interpolate(reshape(S[:,k],origin.shape...), origin, target)
      SInterp[:,k] = vec(A)
    end
    S = SInterp
    grid = target
  elseif !isnothing(fov)
    grid = RegularGridPositions(calibSize(bSF),calibFov(bSF),calibFovCenter(bSF))
  else
    grid = RegularGridPositions(calibSize(bSF),ones(Float64, length(calibSize(bSF))),zeros(Float64, length(calibSize(bSF))))
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

#=
"""
Calculates a noise level from empty measurement `bEmpty`.

    NoiseLevel = getNoiseLevel(bEmpty,bgframes,channels)
"""
function getNoiseLevel(bEmpty, bgframes, channels)
  tmp = getMeasurementsFD(bEmpty, true; frames=bgframes)[:, channels, :, :]
  if length(channels) > 1
    numfreq, numchannels, numpatches, numframes = size(tmp)
    measBG = reshape(permutedims(tmp, [1, 3, 2, 4]), (numfreq * numpatches, numchannels, numframes))
  else
    numchannels = 1
    numfreq, numpatches, numframes = size(tmp)
    measBG = reshape(tmp, (numfreq * numpatches, numchannels, numframes))
  end
  noise = zeros(numfreq * numpatches, numchannels)
  for r = 1:numchannels
    for k = 1:numfreq*numpatches
      tmp = view(measBG, k, r, :)
      #maxBG = mapreduce(abs, max, tmp) # maximum(abs.(measBG[k, r, :]))
      meanBG = mean(tmp)
      noise[k, r] = mean(abs.((tmp .- meanBG)))#./maxBG))
    end
  end
  return mean(noise)
end
=#