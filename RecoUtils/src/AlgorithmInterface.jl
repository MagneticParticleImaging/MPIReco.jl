export AbstractReconstructionAlgorithm, put!, take!, reconstruction

abstract type AbstractReconstructionAlgorithm end

abstract type AbstractReconstructionAlgorithmParameter end

put!(algo::AbstractReconstructionAlgorithm, data) = error("$(typeof(algo)) must implement put!")
take!(algo::AbstractReconstructionAlgorithm) = error("$(typeof(algo)) must implement take!")

function reconstruction(algo::T, u) where {T<:AbstractReconstructionAlgorithm}
  put!(algo, u)
  return take!(algo)
end
