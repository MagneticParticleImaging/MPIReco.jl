export SimilarAlgorithm, SimilarAlgorithmParameter
Base.@kwdef struct SimilarAlgorithmParameter <: AbstractReconstructionAlgorithmParameter
  algo::AbstractReconstructionAlgorithm
end

mutable struct SimilarAlgorithm <: AbstractReconstructionAlgorithm
  params::SimilarAlgorithmParameter
  algo::AbstractReconstructionAlgorithm
  outputChannel::Channel{Any}
end

SimilarAlgorithm(params::SimilarAlgorithmParameter) = SimilarAlgorithm(params, params.algo, Channel{Any}(Inf))

take!(algo::SimilarAlgorithm) = take!(algo.outputChannel)
function put!(algo::SimilarAlgorithm, u)
  result = nothing
  try 
    result = put!(algo.algo, u)
  catch e
    try
      similarAlgo = similar(algo.algo, u)
      result = put!(similarAlgo, u)
      algo.algo = similarAlgo
    catch e2
      throw(e)
    end
  end
  put!(algo.outputChannel, result)
end