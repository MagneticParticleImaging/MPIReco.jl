import Base.length, Base.size

export getSF, SVD, tikhonovLU, setlambda

include("SystemMatrixRecovery.jl")
include("SystemMatrixCenter.jl")
include("SystemMatrixWrapper.jl")

export AbstractSystemMatrixParameter
abstract type AbstractSystemMatrixParameter <: AbstractMPIRecoParameters end

export AbstractSystemMatrixGriddingParameter
abstract type AbstractSystemMatrixGriddingParameter <: AbstractSystemMatrixParameter end
export SystemMatrixGriddingParameter
Base.@kwdef struct SystemMatrixGriddingParameter <: AbstractSystemMatrixGriddingParameter
  gridsize::Vector{Int64} = [1, 1, 1]
  fov::Vector{Float64} = [0.0, 0.0, 0.0]
  center::Vector{Float64} = [0.0,0.0,0.0]
  deadPixels::Vector{Int64} = Int64[]
end
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
end

export SparseSystemMatrixLoadingParameter
Base.@kwdef struct SparseSystemMatrixLoadingParameter{F<:AbstractFrequencyFilterParameter} <: AbstractSystemMatrixLoadingParameter
  freqFilter::F
  sparseTrafo::Union{Nothing, String} = nothing # TODO concrete options here?
  tresh::Float64 = 0.0
  redFactor::Float64 = 0.1
  bgCorrection::Bool = false
  useDFFoV::Bool = false
end
function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, sf::MPIFile, params::DenseSystemMatixLoadingParameter)
  # Construct freqFilter
  freqs = process(t, sf, params.freqFilter)
  S, grid = getSF(sf, freqs, nothing; toKwargs(params)...)
  return freqs, S, grid
end
function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, sf::MPIFile, params::DenseSystemMatixLoadingParameter{<:FrequencyFilteredPreProcessingParameters})
  freqs = params.freqFilter.frequencies
  S, grid = getSF(sf, freqs, nothing; toKwargs(params)...)
  return freqs, S, grid
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
  if solver == "kaczmarz"
    return getSF(bSF, frequencies, sparseTrafo, Kaczmarz; kargs...)
  elseif solver == "pseudoinverse"
    return getSF(bSF, frequencies, sparseTrafo, PseudoInverse; kargs...)
  elseif solver == "cgnr" || solver == "fusedlasso"
    return getSF(bSF, frequencies, sparseTrafo, CGNR; kargs...)
  elseif solver == "direct"
    return getSF(bSF, frequencies, sparseTrafo, DirectSolver; kargs...)
  else
    return getSF(bSF, frequencies, sparseTrafo; kargs...)
  end
end
getSF(bSF, frequencies, sparseTrafo, solver::AbstractLinearSolver; kargs...) = getSF(bSF, frequencies, sparseTrafo, typeof(solver); kargs...)
function getSF(bSF, frequencies, sparseTrafo, solver::Type{<:AbstractLinearSolver}; kargs...)
  SF, grid = getSF(bSF, frequencies, sparseTrafo; kargs...)
  return prepareSF(solver, SF, grid)
end

prepareSF(solver::Type{Kaczmarz}, SF, grid) = transpose(SF), grid
prepareSF(solver::Type{PseudoInverse}, SF, grid) = SVD(svd(transpose(SF))...), grid
prepareSF(solver::Union{Type{CGNR}, Type{FusedLasso}}, SF, grid) = copy(transpose(SF)), grid
prepareSF(solver::Type{DirectSolver}, SF, grid) = RegularizedLeastSquares.tikhonovLU(copy(transpose(SF))), grid
prepareSF(solver::Type{<:RegularizedLeastSquares.AbstractLinearSolver}, SF, grid) = SF, grid



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

function getSF(bSF::MPIFile, frequencies; returnasmatrix = true, procno::Integer=1,
               bgcorrection=false, bgCorrection=bgcorrection, loadasreal=false,
	             gridsize=collect(calibSize(bSF)), numPeriodAverages=1,numPeriodGrouping=1,
	             fov=calibFov(bSF), center=[0.0,0.0,0.0], deadPixels=Int[], kargs...)

  nFreq = rxNumFrequencies(bSF)

  numPeriods = div(acqNumPeriodsPerFrame(bSF),numPeriodGrouping*numPeriodAverages)

  S = getSystemMatrix(bSF, frequencies, bgCorrection=bgCorrection,
                      numPeriodAverages=numPeriodAverages,
                      numPeriodGrouping=numPeriodGrouping; kargs...)

  if !isempty(deadPixels)
    repairDeadPixels(S,gridsize,deadPixels)
  end

  if collect(gridsize) != collect(calibSize(bSF)) ||
    center != [0.0,0.0,0.0] ||
    fov != calibFov(bSF)
    @debug "Perform SF Interpolation..."

    origin = RegularGridPositions(calibSize(bSF),calibFov(bSF),[0.0,0.0,0.0])
    target = RegularGridPositions(gridsize,fov,center)

    SInterp = zeros(eltype(S),prod(gridsize),length(frequencies)*numPeriods)
    for k=1:length(frequencies)*numPeriods
      A = MPIFiles.interpolate(reshape(S[:,k],calibSize(bSF)...), origin, target)
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
