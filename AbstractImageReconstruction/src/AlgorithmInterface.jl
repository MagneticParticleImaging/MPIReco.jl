export AbstractReconstructionAlgorithm
abstract type AbstractReconstructionAlgorithm end

export AbstractReconstructionAlgorithmParameter
abstract type AbstractReconstructionAlgorithmParameter end

export put!, take!
put!(algo::AbstractReconstructionAlgorithm, data) = error("$(typeof(algo)) must implement put!")
take!(algo::AbstractReconstructionAlgorithm) = error("$(typeof(algo)) must implement take!")

export reconstruct
function reconstruct(algo::T, u) where {T<:AbstractReconstructionAlgorithm}
  put!(algo, u)
  return take!(algo)
end

export process
# process(algoT::Type{T}, ...) as pure helper functions
# Overwrite process(algo, ...) to mutate struct based on helper function result
process(algoT::Type{T}, data, param::AbstractReconstructionAlgorithmParameter) where {T<:AbstractReconstructionAlgorithm} = error("No processing defined for algorithm $T with parameter $(typeof(param))")
process(algo::AbstractReconstructionAlgorithm, data, param::AbstractReconstructionAlgorithmParameter) = process(typeof(algo), data, param)

"""
Enable multiple process steps by supplying a Vector of parameters
"""
function process(algo::AbstractReconstructionAlgorithm, data, params::Vector{<:AbstractReconstructionAlgorithmParameter})
  val = process(algo, data, first(params))
  for param ∈ Iterators.drop(params, 1)
    val = process(algo, val, param)
  end
  return val
end

export similar
similar(algo::AbstractReconstructionAlgorithm, data) = error("$(typeof(algo)) must implement similar")
similar(algo::AbstractReconstructionAlgorithm, data, param::AbstractReconstructionAlgorithmParameter) = similar(typeof(algo), data, param)
similar(algoT::Type{T}, data, param::AbstractReconstructionAlgorithmParameter) where {T<:AbstractReconstructionAlgorithm} = deepcopy(param)

export parameter
parameter(algo::AbstractReconstructionAlgorithm) = error("$(typeof(algo)) must implement parameter")