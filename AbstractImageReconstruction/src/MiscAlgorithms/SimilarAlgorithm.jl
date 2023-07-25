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
    put!(algo.algo, u)
    result = take!(algo.algo)
  catch e
    try
      similarAlgo = similar(algo.algo, u)
      put!(similarAlgo, u)
      result = take!(similarAlgo)
      algo.algo = similarAlgo
    catch e2
      throw(e)
    end
  end
  put!(algo.outputChannel, result)
end

parameter(algo::SimilarAlgorithm) = algo.params
similar(algo::SimilarAlgorithm, data) = algo.algo = similar(algo.algo, data)