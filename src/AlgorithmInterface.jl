abstract type AbstractReconstructionAlgorithm end
abstract type AbstractMPIReconstructionAlgorithm <: AbstractReconstructionAlgorithm end

put!(algo::AbstractMPIReconstructionAlgorithm, data) = error("$(typeof(algo)) must implement put!")
take!(algo::AbstractMPIReconstructionAlgorithm) = error("$(typeof(algo)) must implement take!")

function reconstruction(algo::T, u) where T <: AbstractMPIReconstructionAlgorithm
  put!(algo, u)
  return take!(algo)
end

# Traits
abstract type ReconstructionAlgorithmType end

struct SystemMatrixBasedAlgorithm <: ReconstructionAlgorithmType end
struct XSpaceBasedAlgorithm <: ReconstructionAlgorithmType end
struct MachineLearningBasedAlgorithm <: ReconstructionAlgorithmType end
struct MixedAlgorithm <: ReconstructionAlgorithmType end

# Get all child types
# TODO recoAlgorithmType
# TODO undefined for certain "Algorithm" components
recoAlgorithmTypes(::Type{ConcreteRecoAlgorithm}) = SystemMatrixBasedAlgorithm()

# Check if contains
isSystemMatrixBased(::T) where T <: AbstractReconstructionAlgorithm = recoAlgorithmTypes(T) isa SystemMatrixBasedAlgorithm
isXSpaceBased(::T) where T <: AbstractReconstructionAlgorithm = recoAlgorithmTypes(T) isa XSpaceBasedAlgorithm
isMachineLearningBased(::T) where T <: AbstractReconstructionAlgorithm = recoAlgorithmTypes(T) isa MachineLearningBasedAlgorithm
isMixedAlgorithm(::T) where T <: AbstractReconstructionAlgorithm # TODO

abstract type OutputType end
struct IntermediateOutput <: OutputType end
struct ImageOutput <: OutputType

outputType(algo::T) where T <: AbstractMPIReconstructionAlgorithm = IntermediateOutput()
# TODO outputTypes
hasIntermediateOutput(algo::T) where T <: AbstractMPIReconstructionAlgorithm = outputType(algo) isa IntermediateOutput
hasImageOutput(algo::T) where T <: AbstractMPIReconstructionAlgorithm = outputType(algo) isa ImageOutput

include("RuntimeAlgorithms.jl")