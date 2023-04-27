abstract type AbstractMPIReconstructionAlgorithm <: AbstractReconstructionAlgorithm end
abstract type AbstractSystemMatrixReconstructionAlgorithm <: AbstractMPIReconstructionAlgorithm end
abstract type AbstractXSpaceReconstructionAlgorithm <: AbstractMPIReconstructionAlgorithm end
abstract type AbstractMachineLearningReconstructionAlgorithm <: AbstractMPIReconstructionAlgorithm end


abstract type MPIRecoParameters <: AbstractReconstructionAlgorithmParameter end
abstract type AbstractPreProcessingParameters <: MPIRecoParameters end
abstract type AbstractPostProcessingParameters <: MPIRecoParameters end
abstract type AbstractReconstructionParameters <: MPIRecoParameters end
abstract type AbstractRecoAlgorithmParameters <: MPIRecoParameters end
abstract type AbstractBackgroundCorrectionParameters <: MPIRecoParameters end


Base.@kwdef struct MPIFilesPreprocessingParameters <: AbstractPreProcessingParameters
  "Flag whether a background correction shall be applied."
  bgCorrection::Bool = false
  "Flag whether the background frames shall be neglected."
  neglectBGFrames::Bool = true
  "Flag whether a transfer function shall be applied."
  tfCorrection::Bool = true
  "Flag whether the frames shall be sorted (e.g. according to a grid from the positions API)."
  sortFrames::Bool = false
  "Number of frames averages that shall be applied to the raw data. Setting this to `nothing` implies using the default value."
  numAverages::Union{Integer, Nothing} = nothing
  "Number of period averages that shall be applied to the raw data. Setting this to `nothing` implies using the default value."
  numPeriodAverages::Union{Integer, Nothing} = nothing
  "Number of periods that shall be grouped together in the samples dimension. Setting this to `nothing` implies using the default value."
  numPeriodGrouping::Union{Integer, Nothing} = nothing
  # ...
end

# Traits
abstract type ReconstructionAlgorithmType end

struct SystemMatrixBasedAlgorithm <: ReconstructionAlgorithmType end
struct XSpaceBasedAlgorithm <: ReconstructionAlgorithmType end
struct MachineLearningBasedAlgorithm <: ReconstructionAlgorithmType end
struct MixedAlgorithm <: ReconstructionAlgorithmType end

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

include("RuntimeAlgorithms.jl")