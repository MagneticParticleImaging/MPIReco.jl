export AbstractMPIReconstructionAlgorithm
abstract type AbstractMPIReconstructionAlgorithm <: AbstractReconstructionAlgorithm end

export AbstractSystemMatrixReconstructionAlgorithm
abstract type AbstractSystemMatrixReconstructionAlgorithm <: AbstractMPIReconstructionAlgorithm end

export AbstractXSpaceReconstructionAlgorithm
abstract type AbstractXSpaceReconstructionAlgorithm <: AbstractMPIReconstructionAlgorithm end

export AbstractMachineLearningReconstructionAlgorithm
abstract type AbstractMachineLearningReconstructionAlgorithm <: AbstractMPIReconstructionAlgorithm end


export AbstractMPIRecoParameters
abstract type AbstractMPIRecoParameters <: AbstractReconstructionAlgorithmParameter end

export AbstractBackgroundCorrectionParameters
abstract type AbstractBackgroundCorrectionParameters <: AbstractMPIRecoParameters end

export AbstractPreProcessingParameters
abstract type AbstractPreProcessingParameters{T<:AbstractBackgroundCorrectionParameters} <: AbstractMPIRecoParameters end

export AbstractPostProcessingParameters, NoPostProcessing
abstract type AbstractPostProcessingParameters <: AbstractMPIRecoParameters end
struct NoPostProcessing <: AbstractPostProcessingParameters end # TODO remove later
process(algo::AbstractMPIReconstructionAlgorithm, data, ::NoPostProcessing) = data

export AbstractReconstructionParameters
abstract type AbstractReconstructionParameters <: AbstractMPIRecoParameters end

export AbstractRecoAlgorithmParameters
abstract type AbstractRecoAlgorithmParameters <: AbstractMPIRecoParameters end

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
isSystemMatrixBased(::T) where T <: AbstractReconstructionAlgorithm = recoAlgorithmTypes(T) isa SystemMatrixBasedAlgorithm
isXSpaceBased(::T) where T <: AbstractReconstructionAlgorithm = recoAlgorithmTypes(T) isa XSpaceBasedAlgorithm
isMachineLearningBased(::T) where T <: AbstractReconstructionAlgorithm = recoAlgorithmTypes(T) isa MachineLearningBasedAlgorithm
isMixedAlgorithm(::T) where T <: AbstractReconstructionAlgorithm = true # TODO