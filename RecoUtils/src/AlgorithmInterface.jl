export AbstractReconstructionAlgorithm, put!, take!, reconstruction, similar

abstract type AbstractReconstructionAlgorithm end

abstract type AbstractReconstructionAlgorithmParameter end

put!(algo::AbstractReconstructionAlgorithm, data) = error("$(typeof(algo)) must implement put!")
take!(algo::AbstractReconstructionAlgorithm) = error("$(typeof(algo)) must implement take!")

function reconstruction(algo::T, u) where {T<:AbstractReconstructionAlgorithm}
  put!(algo, u)
  return take!(algo)
end

# process(algoT::Type{T}, ...) as pure helper functions
# Overwrite process(algo, ...) to mutate struct based on helper function result
process(algoT::Type{T}, data, param::AbstractReconstructionAlgorithmParameter) where {T<:AbstractReconstructionAlgorithm} = error("No processing defined for algorithm $T with parameter $(typeof(param))")
process(algo::AbstractReconstructionAlgorithm, data, param::AbstractReconstructionAlgorithmParameter) = process(typeof(algo), data, param)

similar(algo::AbstractReconstructionAlgorithm, data) = error("$(typeof(algo)) must implement similar")
similar(algo::AbstractReconstructionAlgorithm, data, param::AbstractReconstructionAlgorithmParameter) = similar(typeof(algo), data, param)
similar(algoT::Type{T}, data, param::AbstractReconstructionAlgorithmParameter) where {T<:AbstractReconstructionAlgorithm} = deepcopy(params)