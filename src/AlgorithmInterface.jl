

abstract type AbstractReconstructionAlgorithm end
abstract type AbstractMPIReconstructionAlgorithm <: AbstractReconstructionAlgorithm end

abstract type AbstractPreprocessingParameters end

Base.@kwdef struct MPIFilesPreprocessingParameters <: AbstractPreprocessingParameters
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

function preprocessingParameters(algo::AbstractReconstructionAlgorithm)
  if :preprocessingParameters in fieldnames(typeof(algo))
    algo.preprocessingParameters
  else
    error("Please implement a method `preprocessingParameters` or add a field `preprocessingParameters` for the algorithm of type `$(typeof(algo))`.")
  end
end

function reconstruction(algo::T, file::MPIFile) where T <: AbstractMPIReconstructionAlgorithm
  if usesStandardPreprocessing(algo)
    params = preprocessingParameters(algo)
    
    splattingDict = Dict{Symbol, Any}()
    for fieldname_ in fieldnames(typeof(params))
      value = getfield(params, fieldname_)
      if !isnothing(value) && fieldname != :neglectBGFrames # `neglectBGFrames` is a non-keyword argument
        splattingDict[fieldname_] = value
      end
    end

    preprocessedData = getMeasurements(file, neglectBGFrames(params); splattingDict...)
  else
    preprocessedData = preprocessing(algo, data)
  end

  reconstructedData = _reconstruction(algo, preprocessedData)

  if hasPostprocessing(algo)
    return postprocessing(algo, reconstructedData)
  else
    return reconstructedData
  end
end

function preprocessing(algo::AbstractReconstructionAlgorithm, file::MPIFile)
  if usesStandardPreprocessing(algo) 
    error("The algorithm uses standard preprocessing from MPIFiles and thus, this method is not implemented for algorithms of type `$(typeof(algo))`.")
  else
    error("The algorithm of type `$(typeof(algo))` does not use standard preprocessing from MPIFiles but the preprocessing method is not defined.")
  end

  return nothing
end

usesStandardPreprocessing(algo::T) where T <: AbstractReconstructionAlgorithm = true
hasPostprocessing(algo::T) where T <: AbstractReconstructionAlgorithm = hasmethod(postprocessing, (T, AbstractArray))

_reconstruction(algo::AbstractReconstructionAlgorithm, preprocessedData::AbstractArray) = error("The algorithm of type `$(typeof(algo))` does not implement a method `_reconstruction`.")

# Traits
abstract type ReconstructionAlgorithmType end

struct SystemMatrixBasedAlgorithm <: ReconstructionAlgorithmType end
struct XSpaceBasedAlgorithm <: ReconstructionAlgorithmType end
struct MachineLearningBasedAlgorithm <: ReconstructionAlgorithmType end

recoAlgorithmType(::Type{ConcreteRecoAlgorithm}) = SystemMatrixBasedAlgorithm()

isSystemMatrixBased(::T) where T <: AbstractReconstructionAlgorithm = recoAlgorithmType(T) isa SystemMatrixBasedAlgorithm
isXSpaceBased(::T) where T <: AbstractReconstructionAlgorithm = recoAlgorithmType(T) isa XSpaceBasedAlgorithm
isMachineLearningBased(::T) where T <: AbstractReconstructionAlgorithm = recoAlgorithmType(T) isa MachineLearningBasedAlgorithm