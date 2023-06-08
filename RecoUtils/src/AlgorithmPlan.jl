export RecoPlan
mutable struct RecoPlan{T<:AbstractReconstructionAlgorithmParameter}
  values::Dict{Symbol, Any}
end

function RecoPlan(::Type{T}; kwargs...) where {T<:AbstractReconstructionAlgorithmParameter}
  kwargs = values(kwargs)
  dict = Dict{Symbol, Any}()
  for field in fieldnames(T)
    dict[field] = haskey(kwargs, field) ? kwargs[field] : missing
  end
  return RecoPlan{T}(dict)
end

Base.propertynames(plan::RecoPlan) = names(plan)

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
Base.setindex!(plan::RecoPlan, name::Symbol, x) = Base.setproperty!(plan, name, x)

Base.ismissing(plan::RecoPlan, name::Symbol) = ismissing(getfield(plan, :values)[name])

export types, type, names
types(::RecoPlan{T}) where {T} = fieldtypes(T)
type(::RecoPlan{T}, name::Symbol) where {T} = fieldtype(T, name)
names(::RecoPlan{T}) where {T} = fieldnames(T)

export build
function build(plan::RecoPlan{T}) where {T}
  nestedPlans = filter(entry -> isa(last(entry), RecoPlan), getfield(plan, :values))
  for (name, nested) in nestedPlans
    plan[name] = build(nested)
  end
  return T(;getfield(plan, :values)...)
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

function toDict!(dict, plan::RecoPlan{T}) where {T}
  dict["type"] = string(T)
  for entry in getfield(plan, :values)
    dict[String(first(entry))] = toDictValue(last(entry))
  end
  return dict
end