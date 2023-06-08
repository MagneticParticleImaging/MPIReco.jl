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
  elseif X <: t || X <: RecoPlan{<:t}
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

toDictType(plan::RecoPlan{T}) where {T} = RecoPlan{getfield(parentmodule(T), nameof(T))}