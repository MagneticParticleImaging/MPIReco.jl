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
export addRecoPlanPath, getRecoPlanList
const recoPlanPaths = [normpath(joinpath(@__DIR__, "..", "config"))]
"""
    addRecoPlanPath(path::String)

Add all `RecoPlans` within the given directory as potential 
"""
addRecoPlanPath(path::String) = !(path in recoPlanPaths) ? pushfirst!(recoPlanPaths, path) : nothing
"""
    getRecoPlanList(; full = false)

Retrieve a list of currently known `RecoPlan`s. If `full` is `true` then full paths are returned. 
"""
function getRecoPlanList(; full = false)
  result = String[]
  for path in recoPlanPaths
    if isdir(path)
      push!(result, [plan for plan in filter(a->contains(a,".toml"),readdir(path, join =  full))]...)
    end
  end
  return result
end



function planpath(name::AbstractString)
  if isfile(name)
    return name
  end

  for dir in recoPlanPaths
    filename = joinpath(dir, string(name, ".toml"))
    if isfile(filename)
      return filename
    end
  end
  throw(ArgumentError("Could not find a suitable MPI reconstruction plan with name $name. Custom plans can be stored in the following directories $(join(recoPlanPaths, ", "))."))
end

const recoPlans = LRU{UInt64, RecoPlan}(maxsize = 3)

export reconstruct
"""
    reconstruct(name::AbstractString, data::MPIFile, cache::Bool = false, modules = [AbstractImageReconstruction, MPIFiles, MPIReco, RegularizedLeastSquares]; kwargs...)

Perform a reconstruction with the `RecoPlan` specified by `name` and given `data`.
Additional keyword arguments can be passed to the reconstruction plan.

`RecoPlans` can be stored in the in the MPIReco package config folder or in a custom folder. New folder can be added with `addRecoPlanPath`. The first plan found is used.
Alternatively, name can be a path to specific plan file.

If `cache` is `true` the reconstruction plan is cached and reused if the plan file has not changed. If a keyword argument changes the structure of the plan the cache is also bypassed.
The cache considers the last modification time of the plan file and can be manually be emptied with `emptyRecoCache!()`.

# Examples
```julia
julia> mdf = MPIFile("data.mdf");

julia> reconstruct("SinglePatch", mdf; solver = Kaczmarz, reg = [L2Regularization(0.3f0)], iterations = 10, frames = 1:10, ...)
```
"""
function reconstruct(name::AbstractString, data::MPIFile, cache::Bool = false, modules = [AbstractImageReconstruction, MPIFiles, MPIReco, RegularizedLeastSquares]; kwargs...)
  plan = loadRecoPlan(name, cache, modules; kwargs...)
  setAll!(plan; kwargs...)
  return reconstruct(build(plan), data)
end
function loadRecoPlan(name::AbstractString, cache::Bool, modules; kwargs...)
  planfile = planpath(name)

  # If the user disables caching or changes the plan structure we bypass the cache
  kwargValues = values(values(kwargs))
  if !cache || any(val -> isa(val, AbstractRecoPlan) || isa(val, AbstractImageReconstructionParameters), kwargValues)
    return loadPlan(planfile, modules)
  end

  key = hash(planfile, hash(mtime(planfile)))
  return get!(recoPlans, key) do
    loadPlan(planfile, modules)
  end
end

export emptyRecoCache!
"""
    emptyRecoCache!()

Empty the cache of `RecoPlans`. This is useful if the cache is too large.
"""
emptyRecoCache!() = Base.empty!(recoPlans)

# Check if contains
isSystemMatrixBased(::T) where T <: AbstractImageReconstructionAlgorithm = recoAlgorithmTypes(T) isa SystemMatrixBasedAlgorithm
isXSpaceBased(::T) where T <: AbstractImageReconstructionAlgorithm = recoAlgorithmTypes(T) isa XSpaceBasedAlgorithm
isMachineLearningBased(::T) where T <: AbstractImageReconstructionAlgorithm = recoAlgorithmTypes(T) isa MachineLearningBasedAlgorithm
isMixedAlgorithm(::T) where T <: AbstractImageReconstructionAlgorithm = true # TODO