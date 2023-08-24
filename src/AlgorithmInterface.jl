export AbstractMPIRecoAlgorithm
abstract type AbstractMPIRecoAlgorithm <: AbstractImageReconstructionAlgorithm end

export AbstractSystemMatrixMPIRecoAlgorithm
abstract type AbstractSystemMatrixMPIRecoAlgorithm <: AbstractMPIRecoAlgorithm end

export AbstractXSpaceMPIRecoAlgorithm
abstract type AbstractXSpaceMPIRecoAlgorithm <: AbstractMPIRecoAlgorithm end

export AbstractMachineLearningMPIRecoAlgorithm
abstract type AbstractMachineLearningMPIRecoAlgorithm <: AbstractMPIRecoAlgorithm end


export AbstractMPIRecoParameters
abstract type AbstractMPIRecoParameters <: AbstractImageReconstructionParameters end

export AbstractMPIBackgroundCorrectionParameters
abstract type AbstractMPIBackgroundCorrectionParameters <: AbstractMPIRecoParameters end

export AbstractMPIPreProcessingParameters
abstract type AbstractMPIPreProcessingParameters{T<:AbstractMPIBackgroundCorrectionParameters} <: AbstractMPIRecoParameters end

export AbstractMPIPostProcessingParameters, NoPostProcessing
abstract type AbstractMPIPostProcessingParameters <: AbstractMPIRecoParameters end
struct NoPostProcessing <: AbstractMPIPostProcessingParameters end # TODO remove later
process(algo::AbstractMPIRecoAlgorithm, data, ::NoPostProcessing) = data

export AbstractMPIReconstructionParameters
abstract type AbstractMPIReconstructionParameters <: AbstractMPIRecoParameters end

export AbstractMPIRecoAlgorithmParameters
abstract type AbstractMPIRecoAlgorithmParameters <: AbstractMPIRecoParameters end

# Traits
abstract type ReconstructionAlgorithmType end

struct SystemMatrixBasedAlgorithm <: ReconstructionAlgorithmType end
struct XSpaceBasedAlgorithm <: ReconstructionAlgorithmType end
struct MachineLearningBasedAlgorithm <: ReconstructionAlgorithmType end
struct MixedAlgorithm <: ReconstructionAlgorithmType end

# TODO recoAlgorithmType
# TODO undefined for certain "Algorithm" components
#recoAlgorithmTypes(::Type{ConcreteRecoAlgorithm}) = SystemMatrixBasedAlgorithm()
export plandir
plandir() = abspath(homedir(), ".mpi", "RecoPlans")

# Check if contains
isSystemMatrixBased(::T) where T <: AbstractImageReconstructionAlgorithm = recoAlgorithmTypes(T) isa SystemMatrixBasedAlgorithm
isXSpaceBased(::T) where T <: AbstractImageReconstructionAlgorithm = recoAlgorithmTypes(T) isa XSpaceBasedAlgorithm
isMachineLearningBased(::T) where T <: AbstractImageReconstructionAlgorithm = recoAlgorithmTypes(T) isa MachineLearningBasedAlgorithm
isMixedAlgorithm(::T) where T <: AbstractImageReconstructionAlgorithm = true # TODO