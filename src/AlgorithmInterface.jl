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
process(algo::AbstractMPIRecoAlgorithm, ::NoPostProcessing, data) = data

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
export planpath, plandir
plandir() = abspath(homedir(), ".mpi", "RecoPlans")

function planpath(name::AbstractString)
  for dir in [joinpath(@__DIR__, "..", "config"),plandir()]
    filename = joinpath(dir, string(name, ".toml"))
    if isfile(filename)
      return filename
    end
  end
  throw(ArgumentError("Could not find a suitable MPI reconstruction plan with name $name.\nCustom plans can be stored in $(plandir())."))
end

const recoPlans = LRU{UInt64, RecoPlan}(maxsize = 3)

export reconstruct
function reconstruct(name::AbstractString, data::MPIFile, cache::Bool = true; kwargs...) 
  planfile = AbstractImageReconstruction.planpath(MPIReco, name)
  key = hash(planfile, hash(mtime(planfile)))
  plan = get!(recoPlans, key) do
    loadPlan(planfile, [AbstractImageReconstruction, MPIFiles, MPIReco, RegularizedLeastSquares])
  end
  setAll!(plan; kwargs...)
  return reconstruct(build(plan), data)
end

export emptyRecoCache!
emptyRecoCache!() = Base.empty!(recoPlans)

# Check if contains
isSystemMatrixBased(::T) where T <: AbstractImageReconstructionAlgorithm = recoAlgorithmTypes(T) isa SystemMatrixBasedAlgorithm
isXSpaceBased(::T) where T <: AbstractImageReconstructionAlgorithm = recoAlgorithmTypes(T) isa XSpaceBasedAlgorithm
isMachineLearningBased(::T) where T <: AbstractImageReconstructionAlgorithm = recoAlgorithmTypes(T) isa MachineLearningBasedAlgorithm
isMixedAlgorithm(::T) where T <: AbstractImageReconstructionAlgorithm = true # TODO