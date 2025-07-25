export LowLevelReconstructionAlgorithm
Base.@kwdef struct LowLevelReconstructionAlgorithm{P <: LeastSquaresParameters} <: AbstractMPIRecoAlgorithm
  params::P
  output::Channel{Any}
  LowLevelReconstructionAlgorithm(params::P) where P = new{P}(params, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{LowLevelReconstructionAlgorithm}) = SystemMatrixBasedAlgorithm()
AbstractImageReconstruction.parameter(algo::LowLevelReconstructionAlgorithm) = algo.params
Base.lock(algo::LowLevelReconstructionAlgorithm) = lock(algo.output)
Base.unlock(algo::LowLevelReconstructionAlgorithm) = unlock(algo.output)
Base.isready(algo::LowLevelReconstructionAlgorithm) = isready(algo.output)
Base.wait(algo::LowLevelReconstructionAlgorithm) = wait(algo.output)
AbstractImageReconstruction.take!(algo::LowLevelReconstructionAlgorithm) = Base.take!(algo.output)

function AbstractImageReconstruction.put!(algo::LowLevelReconstructionAlgorithm, u::AbstractArray)
  put!(algo.output, process(algo, algo.params, u))
end
