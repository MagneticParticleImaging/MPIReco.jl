# Can be shared across algorithms
Base.@kwdef struct CommonPreProcessingParameters <: AbstractPreProcessingParameter
  bgCorrection::Bool = false
  numAverages::Int64 = 1
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  recChannels::UnitRange{Int64} = 1:rxNumChannels(sf)
end

# These one can be more tailored to one algo
# Or common system matrix parameter
Base.@kwdef struct SinglePatchReconstructionParameter{S<:AbstractLinearSolver} <: AbstractReconstructionParameters
  # File
  sf::MPIFile
  loadasreal::Bool = false
  sparseTrafo::Union{Nothing, Any} = nothing
  #saveTrafo::Bool = false
  # Freqs
  minFreq::Int64 = 0
  maxFreq::Int64 = rxBandwidth(sf)
  numUsedFreqs::Int64 = -1
  # SNR
  threshSNR::Float64=-1
  sortBySNR::Bool = false
  # Solver
  solver::Type{S}
  # Grid
  gridsize::Vector{Int64} = gridSizeCommon(bSF)
  fov::Vector{Float64} = calibFov(bSF)
  center::Vector{Float64} = [0.0,0.0,0.0]
  useDFFoV::Bool = false
  # Misc.
  frames::Union{Nothing, UnitRange{Int64}} = nothing
  deadPixels::Vector{Int64} = Int64[]
end

# Having the parameters for an algorithm in nested structs allows for easier introspection
# Using AbstractParameters here allows one to combine different parameters to potentially achieve different algorithms
# See more in later comments.
# Could also make this a parametric type, though it shouldnt be performance critical
Base.@kwdef mutable struct SinglePatchParameters <: AbstractRecoAlgorithmParameters
  pre::Union{AbstractPreProcessingParameters,Nothing} = nothing
  reco::Union{AbstractReconstructionParameters,Nothing} = nothing
  post::Union{AbstractPostProcessingParameters,Nothing} = nothing
end

Base.@kwdef mutable struct SinglePatchReconstruction <: AbstractMPIReconstructionAlgorithm
  params::SinglePatchParameters
  # Could also do reconstruction progress meter here
  S
  grid
  freqs
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchReconstructionParameter)
  # Prepare system matrix based on pre and reco params
  freqs = getFreqs(params.pre, params.reco)
  S, grid = getSF(params.pre, params.reco, freqs)
  
  # Also possible depending on the reco params to adapt and change the pre/post params
  # User might input common pre params, but if they want temp. reg., we could build specifc
  # param structs. Or param structs could be user friendly such as min/max freq, and here we compute the proper freqs array

  # It is also possible to nest reconstruction algorithms here

  # Something that is still missing is, that an algorithm needs to be able to encode
  # which parameter combinations are possible or this needs to be checked here
  # Similar to requiredDevices in MPIMeasurements.

  # At the end we want to be in a state where we can handle a series of put! requests
  return SinglePatchReconstruction(params, S, grid, freqs, Channel{Any}(32))
end

getFreqs(pre::CommonPreProcessingParameters, reco::SinglePatchReconstructionParameter) = filterFrequences(reco.sf, toKargs([pre, reco])...)

# TODO Maybe do this on Type{T}
function getSF(pre::CommonPreProcessingParameters, reco::SinglePatchReconstructionParameter, freqs)
  return getSF(reco.sf, freqs, reco.sparseTrafo, reco.solver; toKwargs([pre, reco])...)
end

take!(algo::SinglePatchReconstruction) = take!(algo.output)

# Similar to the already propossed generic reconstruction, but that's just this example
# We could for example have a reco where reconstruction + post loops until a condition is met
# This way we dont need to map these things on generic pre, reco or post steps/methods

# However with this more "free-form" put! method, we might have more code repetition
# The idea here is to provide a lot of flexible and reusable "helper" functions with dispatch on parameters
function put!(algo::SinglePatchReconstruction, data)
  consistenceCheck(bSF, bMeas)
  pre = process(algo, data, algo.params.pre)
  result = process(algo, pre, algo.params.reco)
  if !isnothing(algo.params.post)
    result = process(algo, result, algo.params.post)
  end
  put!(algo.output, result)
end

function similar(algo::SinglePatchReconstruction, data::MPIFile)
  pre = similar(algo, data, algo.param.pre)
  reco = similar(algo, data, algo.param.reco)
  post = similar(algo, data, algo.parms.post)

  # This isnt a nice structure atm, because the corrections of the params cant be done individually
  # and in certain parts should only be applied after construction

  # Check if freq still valid
  freqs = algo.freqs
  S = algo.S
  grid = algo.grid
  if rxNumSamplingPoints(reco.sf) == rxNumSamplingPoints(data)
    # Ensure that no frequencies are used that are not present in the measurement
    freqs = intersect(algo.freqs, filterFrequencies(data, toKwargs(pre)))
    if freqs != algo.freqs
      S, grid = getSF(params.pre, params.reco, freqs)
    end
  end

  result = SinglePatchReconstruction(SinglePatchReconstructionParameter(pre, reco, post), S, grid, freqs, Channel{Any}(32))

  # TODO this only works for common preprocessing, because i directly access fields
  # If S is processed and fits not to the measurements because of numPeriodsGrouping
  # or numPeriodAverages being applied we need to set these so that the 
  # measurements are loaded correctly
  if rxNumSamplingPoints(bSF) > rxNumSamplingPoints(bMeas)
    pre.numPeriodGrouping = rxNumSamplingPoints(bSF) ÷ rxNumSamplingPoints(bMeas)
  end
  if acqNumPeriodsPerFrame(bSF) < acqNumPeriodsPerFrame(bMeas)
    pre.numPeriodAverages = acqNumPeriodsPerFrame(bMeas) ÷ (acqNumPeriodsPerFrame(bSF) * numPeriodGrouping)
  end
  return result
end

# These dont need to be specified for each algorithm, but each algorithm can have unique behaviour
function process(algo::AbstractMPIReconstructionAlgorithm, f::MPIFile, params::CommonPreProcessingParameters)
  result = getMeasurementFD(f, ...) 
end

function process(algo::AbstractMPIReconstructionAlgorithm, u::Array, params::CommonPreProcessingParameters)
  result = ...
end

# Alternatively one might do something like
foo(algo::Type{AbstractMPIReconstructionAlgorithm}, ...)
# If one wants to avoid generic functions interacting with the algo struct
# Type is just used like a trait for dispatch