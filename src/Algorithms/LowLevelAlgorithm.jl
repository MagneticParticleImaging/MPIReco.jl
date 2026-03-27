export LowLevelReconstructionAlgorithm
@reconstruction struct LowLevelReconstructionAlgorithm{P <: LeastSquaresParameters} <: AbstractMPIRecoAlgorithm
  @parameter params::P
end
recoAlgorithmTypes(::Type{LowLevelReconstructionAlgorithm}) = SystemMatrixBasedAlgorithm()