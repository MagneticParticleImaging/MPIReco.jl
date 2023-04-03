# Can be shared across algorithms
Base.@kwdef struct CommonPreProcessingParameters <: AbstractPreProcessingParameter
  # ...
  "Flag whether a background correction shall be applied."
  bgCorrection::Bool = false
  "Number of frames averages that shall be applied to the raw data. Setting this to `nothing` implies using the default value."
  numAverages::Union{Integer,Nothing} = nothing
  "Number of period averages that shall be applied to the raw data. Setting this to `nothing` implies using the default value."
  numPeriodAverages::Union{Integer,Nothing} = nothing
  "Number of periods that shall be grouped together in the samples dimension. Setting this to `nothing` implies using the default value."
  numPeriodGrouping::Union{Integer,Nothing} = nothing
  # ...
end

Base.@kwdef struct CommonPostProcessingParameters <: AbstractPostProcessingParameter
  # ...
end

# These one can be more tailored to one algo
Base.@kwdef struct SinglePatchReconstructionParameter <: AbstractReconstructionParameters
  # ...
  # SystemMatrix
  # tempRegularization::Bool = false
  # solver
  # lambda
  # ...
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
  # S, grid = ...
  # Compute freqs 
  # freqs = 
  
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

take!(algo::SinglePatchReconstruction) = take!(algo.output)

# Similar to the already propossed generic reconstruction, but that's just this example
# We could for example have a reco where reconstruction + post loops until a condition is met
# This way we dont need to map these things on generic pre, reco or post steps/methods

# However with this more "free-form" put! method, we might have more code repetition
# The idea here is to provide a lot of flexible and reusable "helper" functions with dispatch on parameters
function put!(algo::SinglePatchReconstruction, data)
  pre = process(algo, f, algo.params.pre)
  result = reconstruct(algo, pre, algo.params.reco)
  if !isnothing(algo.params.post)
    result = process(algo, result, algo.params.post)
  end
  put!(algo.output, result)
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