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
const DEFAULT_PLANS_PATH = @path joinpath(@__DIR__, "..", "config")
const recoPlanPaths = AbstractString[DEFAULT_PLANS_PATH]
"""
    addRecoPlanPath(path)

Add all `RecoPlans` within the given directory as potential 
"""
addRecoPlanPath(path) = !(path in recoPlanPaths) ? pushfirst!(recoPlanPaths, path) : nothing
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

export addRecoPlanModule, getRecoPlanModules
const recoPlanModules::Vector{Module} = [AbstractImageReconstruction, MPIFiles, MPIReco, RegularizedLeastSquares]
"""
    addRecoPlanPath(mod::Module)

Add the `mod` module to the list of modules which are used for plan loading 
"""
addRecoPlanModule(mod::Module) = !(mod in recoPlanModules) ? push!(recoPlanModules, mod) : nothing
"""
    getRecoPlanModules()

Retrieve a list of currently used modules for plan loading 
"""
getRecoPlanModules() = recoPlanModules

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
    reconstruct(name::AbstractString, data::MPIFile, cache::Bool = false, modules = getRecoPlanModules(); kwargs...)

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
function reconstruct(name::AbstractString, data::MPIFile, cache::Bool = false, modules = getRecoPlanModules(); kwargs...)
  plan = loadRecoPlan(name, cache, modules; kwargs...)
  return reconstruct(build(plan), data)
end
# Load plan with RecoCache consideration
function loadRecoPlan(name::AbstractString, cache::Bool, modules; kwargs...)
  planfile = planpath(name)

  # If the user disables caching or changes the plan structure we bypass the cache
  kwargValues = values(values(kwargs))
  if !cache || any(val -> isa(val, AbstractRecoPlan) || isa(val, AbstractImageReconstructionParameters), kwargValues)
    plan = loadRecoPlan(planfile, modules; kwargs...)
  else 
    key = hash(planfile, hash(mtime(planfile)))
    plan = get!(recoPlans, key) do
      loadRecoPlan(planfile, modules; kwargs...)
    end
  end

  return plan
end
# Load plan from a .toml file
function loadRecoPlan(planfile::AbstractString, modules; kwargs...)
  return open(planfile, "r") do io
    return loadRecoPlan(io, modules; kwargs...)
  end
end
# Load plan from an io (could be file or iobuffer backed string)
function loadRecoPlan(io, modules; kwargs...)
  plan = loadPlan(io, modules)
  setFirst = filter(kw->kw[2] isa Union{<:AbstractImageReconstructionAlgorithm, <:AbstractImageReconstructionParameters, <:AbstractRecoPlan}, kwargs)
  setAll!(plan; setFirst...)
  setAll!(plan; [kw for kw in kwargs if kw ∉ setFirst]...)
  return plan
end
export MPIRecoPlan
"""
    MPIRecoPlan(value, modules = getRecoPlanModules(); kwargs...)

Load a `RecoPlan` and set the keyword arguments if applicable. Value can be the name of a plan or path to plan file or an MDF or a path to one.
"""
function MPIRecoPlan(value::AbstractString, modules = getRecoPlanModules(); kwargs...)
  if isfile(value) && endswith(value, ".toml")
    return loadRecoPlan(value, modules; kwargs...)
  elseif isfile(value)
    return MPIRecoPlan(MDFFile(value), modules; kwargs...)
  else
    return loadRecoPlan(planpath(value), modules; kwargs...)
  end
end
function MPIRecoPlan(file::MDFFileV2, modules; kwargs...)
  error("Loading reconstruction from MDF not yet supported")
  planstr = ""
  io = IOBuffer(planstr)
  return loadRecoPlan(io, modules; kwargs...) 
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