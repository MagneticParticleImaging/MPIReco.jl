export RecoPlan
mutable struct RecoPlan{T<:Union{AbstractReconstructionAlgorithmParameter, AbstractReconstructionAlgorithm}}
  values::Dict{Symbol, Any}
end

function RecoPlan(::Type{T}; kwargs...) where {T<:AbstractReconstructionAlgorithmParameter}
  kwargs = values(kwargs)
  dict = Dict{Symbol, Any}()
  for field in fieldnames(T)
    dict[field] = haskey(kwargs, field) ? kwargs[field] : missing
  end
  return RecoPlan{getfield(parentmodule(T), nameof(T))}(dict)
end
function RecoPlan(::Type{T}; kwargs...) where {T<:AbstractReconstructionAlgorithm}
  dict = Dict{Symbol, Any}()
  mod = parentmodule(T)
  dict[:module] = mod
  dict[:version] = string(pkgversion(mod))
  dict[:parameter] = missing
  return RecoPlan{getfield(parentmodule(T), nameof(T))}(dict)
end



Base.propertynames(plan::RecoPlan{T}) where {T} = keys(getfield(plan, :values))

Base.getproperty(plan::RecoPlan{T}, name::Symbol) where {T} = getfield(plan, :values)[name]
Base.getindex(plan::RecoPlan{T}, name::Symbol) where {T} = Base.getproperty(plan, name)

function Base.setproperty!(plan::RecoPlan{T}, name::Symbol, x::X) where {T, X}
  t = type(plan, name)
  if !haskey(getfield(plan, :values), name)
    error("type $T has no field $name")
  elseif X <: t || X <: RecoPlan{<:t} || ismissing(x)
    getfield(plan, :values)[name] = x
  else
    getfield(plan, :values)[name] = convert(t, x)
  end
  return Base.getproperty(plan, name)
end
Base.setindex!(plan::RecoPlan, x, name::Symbol) = Base.setproperty!(plan, name, x)

Base.ismissing(plan::RecoPlan, name::Symbol) = ismissing(getfield(plan, :values)[name])

export types, type
types(::RecoPlan{T}) where {T<:AbstractReconstructionAlgorithmParameter} = fieldtypes(T)
type(::RecoPlan{T}, name::Symbol) where {T<:AbstractReconstructionAlgorithmParameter} = fieldtype(T, name)

function type(plan::RecoPlan{T}, name::Symbol) where {T<:AbstractReconstructionAlgorithm}
  if name == :version
    return String
  elseif name == :parameter
    return Union{Missing, RecoPlan}
  elseif name == :module
    return Module
  else
    error("type $(typeof(plan)) has no field $name")
  end
end
types(::RecoPlan{T}) where {T<:AbstractReconstructionAlgorithm} = [type(plan, name) for name in propertynames(plan)]


export clear!
function clear!(plan::RecoPlan{T}, preserve::Bool = true) where {T<:AbstractReconstructionAlgorithmParameter}
  dict = getfield(plan, :values)
  for key in keys(dict)
    value = dict[key]
    dict[key] = typeof(value) <: RecoPlan && preserve ? clear!(value, preserve) : missing
  end
  return plan
end
clear!(plan::RecoPlan{T}, preserve::Bool = true) where {T<:AbstractReconstructionAlgorithm} = clear!(plan.parameter, preserve)

export build
function build(plan::RecoPlan{T}) where {T<:AbstractReconstructionAlgorithmParameter}
  fields = getfield(plan, :values)
  nestedPlans = filter(entry -> isa(last(entry), RecoPlan), fields)
  for (name, nested) in nestedPlans
    plan[name] = build(nested)
  end
  fields = filter(entry -> !ismissing(last(entry)), fields)
  return T(;fields...)
end
function build(plan::RecoPlan{T}) where {T<:AbstractReconstructionAlgorithm}
  parameter = build(plan[:parameter])
  return T(parameter)
end

export toPlan
function toPlan(param::AbstractReconstructionAlgorithmParameter)
  args = Dict{Symbol, Any}()
  for field in fieldnames(typeof(param))
    value = getproperty(param, field)
    if typeof(value) <: AbstractReconstructionAlgorithmParameter
      args[field] = toPlan(value)
    else
      args[field] = value
    end
  end
  return RecoPlan(typeof(param); args...)
end
toPlan(algo::AbstractReconstructionAlgorithm, params::AbstractReconstructionAlgorithmParameter) = toPlan(typeof(algo), params) 
function toPlan(::Type{T}, params::AbstractReconstructionAlgorithmParameter) where {T<:AbstractReconstructionAlgorithm}
  plan = RecoPlan(T)
  plan[:parameter] = toPlan(params)
  return plan
end

toDictModule(plan::RecoPlan{T}) where {T} = parentmodule(T)
toDictType(plan::RecoPlan{T}) where {T} = RecoPlan{getfield(parentmodule(T), nameof(T))}
function toDictValue!(dict, value::RecoPlan, field::Symbol)
  x = getproperty(value, field)
  if !ismissing(x)
    dict[string(field)] = toDictValue(x)
  end
  return dict
end

export loadPlan
function loadPlan(filename::AbstractString, modules::Vector{Module})
  dict = TOML.parsefile(filename)
  modDict = createModuleDataTypeDict(modules)
  return loadPlan!(dict, modDict)
end
function createModuleDataTypeDict(modules::Vector{Module})
  modDict = Dict{String, Dict{String, Union{DataType, UnionAll}}}()
  for mod in modules
    typeDict = Dict{String, Union{DataType, UnionAll}}()
    for field in names(mod)
      try
        t = getfield(mod, field)
        if t isa DataType || t isa UnionAll
          typeDict[string(t)] = t
        end
      catch
      end
    end
    modDict[string(mod)] = typeDict
  end
  return modDict
end
function loadPlan!(dict::Dict{String, Any}, modDict::Dict{String, Dict{String, Union{DataType, UnionAll}}})
  re = r"RecoPlan\{(.*)\}"
  m = match(re, pop!(dict, ".type"))
  if !isnothing(m)
    type = m.captures[1]
    mod = pop!(dict, ".module")
    plan = RecoPlan(modDict[mod][type])
    loadPlan!(plan, dict, modDict)
    return plan
  else
    # Has to be parameter or algo or broken toml
  end
end
function loadPlan!(plan::RecoPlan{T}, dict::Dict{String, Any}, modDict::Dict{String, Dict{String, Union{DataType, UnionAll}}}) where {T<:AbstractReconstructionAlgorithmParameter}
  for name in propertynames(plan)
    t = type(plan, name)
    param = missing
    key = string(name)
    if haskey(dict, key)
      if t <: AbstractReconstructionAlgorithm || t <: AbstractReconstructionAlgorithmParameter
        param = loadPlan!(dict[key], modDict)
      # Handle Type{T}, structs with 0 fields and structs with more (assume kwargs for later?)
      # I think for some of these we need to go back into the modDict to further specify a type
      # As the Plan type can be abstract
      else
        param = loadPlanValue!(t, dict[key], modDict)
      end
    end
    plan[name] = param
  end
end
# Type{<:T} where {T} fields
function loadPlanValue!(::UnionAll, value::Dict, modDict::Dict{String, Dict{String, Union{DataType, UnionAll}}})
  re = r"Type\{(.*)\}"
  m = match(re, pop!(value, ".type"))
  mod = pop!(value, ".module")
  type = m.captures[1]
  return modDict[mod][type]
end
# Vectors
function loadPlanValue!(::Type{Vector{<:T}}, value::Vector, modDict) where {T}
  result = Any[]
  for val in value
    type = modDict[pop!(val, ".module")][pop!(val, ".type")]
    push!(result, fromTOML(type, val))
  end
  # Narrow vector
  return identity.(result)
end
# Structs with fields
function loadPlanValue!(t::Type{T}, value::Dict, modDict) where {T}
  if !in(value[".module"], keys(modDict))
    return fromTOML(t, value)
  else
    type = modDict[pop!(value, ".module")][pop!(value, ".type")]
    kwargs = Dict(Symbol(x) => dict[x] for x in keys(value))
    return type(kwargs...)
  end
end
loadPlanValue!(t, value, modDict) = fromTOML(t, value)

export savePlan
savePlan(filename::AbstractString, plan::RecoPlan) = toTOML(filename, plan)