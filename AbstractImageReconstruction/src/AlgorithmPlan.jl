export RecoPlan

export AbstractPlanListener, TransientListener, SerializableListener
abstract type AbstractPlanListener end
abstract type TransientListener <: AbstractPlanListener end
abstract type SerializableListener <: AbstractPlanListener end

mutable struct RecoPlan{T<:Union{AbstractReconstructionAlgorithmParameter, AbstractReconstructionAlgorithm}}
  parent::Union{Nothing, RecoPlan}
  values::Dict{Symbol, Any}
  listeners::Dict{Symbol, Vector{AbstractPlanListener}}
  setProperties::Dict{Symbol, Bool}
  function RecoPlan(::Type{T}; kwargs...) where {T<:AbstractReconstructionAlgorithmParameter}
    dict = Dict{Symbol, Any}()
    listeners = Dict{Symbol, Vector{AbstractPlanListener}}()
    setProperties = Dict{Symbol, Bool}()
    for field in fieldnames(T)
      dict[field] =  missing
      listeners[field] = AbstractPlanListener[]
      setProperties[field] = false
    end
    plan = new{getfield(parentmodule(T), nameof(T))}(nothing, dict, listeners, setProperties)
    setvalues!(plan, kwargs...)
    return plan
  end
  function RecoPlan(::Type{T}) where {T<:AbstractReconstructionAlgorithm}
    dict = Dict{Symbol, Any}()
    listeners = Dict{Symbol, Vector{AbstractPlanListener}}()
    setProperties = Dict{Symbol, Bool}()
    dict[:parameter] = missing
    listeners[:parameter] = AbstractPlanListener[]
    setProperties[:parameter] = false
    return new{getfield(parentmodule(T), nameof(T))}(nothing, dict, listeners, setProperties)
  end  
end

function setvalues!(plan::RecoPlan{T}; kwargs...) where {T<:AbstractReconstructionAlgorithmParameter}
  kwargs = values(kwargs)
  for field in propertynames(plan)
    if haskey(kwargs, field)
      setvalue!(plan, field, kwargs[field])
    end
  end
end

export parent, parent!
parent(plan::RecoPlan) = getfield(plan, :parent)
parent!(plan::RecoPlan, parent::RecoPlan) = setfield!(plan, :parent, parent)
function parentfields(plan::RecoPlan)
  trace = Symbol[]
  return parentfields!(trace, plan)
end
function parentfields!(trace::Vector{Symbol}, plan::RecoPlan)
  p = parent(plan)
  if !isnothing(p)
    for property in propertynames(p)
      if getproperty(p, property) === plan
        pushfirst!(trace, property)
        return parentfields!(trace, p)
      end
    end
  end
  return trace
end

function RecoPlan(parent::RecoPlan, t::Type; kwargs...)
  plan = RecoPlan(t; kwargs...)
  parent!(plan, parent)
  return plan
end

Base.propertynames(plan::RecoPlan{T}) where {T} = keys(getfield(plan, :values))

Base.getproperty(plan::RecoPlan{T}, name::Symbol) where {T} = getfield(plan, :values)[name]
Base.getindex(plan::RecoPlan{T}, name::Symbol) where {T} = Base.getproperty(plan, name)

export propertyupdate!, ispropertyset, setvalue!, valueupdate
function propertyupdate!(listener::AbstractPlanListener, origin, field, old, new)
  # NOP
end
function Base.setproperty!(plan::RecoPlan{T}, name::Symbol, x::X) where {T, X}
  old = getproperty(plan, name) 
  setvalue!(plan, name, x)
  getfield(plan, :setProperties)[name] = true
  for listener in getlisteners(plan, name)
    try
      propertyupdate!(listener, plan, name, old, x)
    catch e
      @error "Exception in listener $listener " e
    end
  end
end
ispropertyset(plan::RecoPlan, name::Symbol) = getfield(plan, :setProperties)[name]
Base.setindex!(plan::RecoPlan, x, name::Symbol) = Base.setproperty!(plan, name, x)
function valueupdate(listener::AbstractPlanListener, origin, field, old, new)
  # NOP
end
function setvalue!(plan::RecoPlan{T}, name::Symbol, x::X) where {T, X}
  old = Base.getproperty(plan, name)
  t = type(plan, name)
  if !haskey(getfield(plan, :values), name)
    error("type $T has no field $name")
  elseif X <: t || X <: RecoPlan{<:t} || ismissing(x)
    getfield(plan, :values)[name] = x
  else
    getfield(plan, :values)[name] = convert(t, x)
  end
  new = Base.getproperty(plan, name)
  for listener in getlisteners(plan, name)
    try
      valueupdate(listener, plan, name, old, new)
    catch e
      @error "Exception in listener $listener " e
    end
  end
  return new
end

Base.ismissing(plan::RecoPlan, name::Symbol) = ismissing(getfield(plan, :values)[name])

export setAll!
function setAll!(plan::RecoPlan{T}, name::Symbol, x) where {T<:AbstractReconstructionAlgorithmParameter}
  fields = getfield(plan, :values)
  nestedPlans = filter(entry -> isa(last(entry), RecoPlan), fields)
  for (key, nested) in nestedPlans
    key != name && setAll!(nested, name, x)
  end
  if haskey(fields, name)
    try
      plan[name] = x
    catch ex
      @warn "Could not set $name of $T with value of type $(typeof(x))"
    end
  end
end
setAll!(plan::RecoPlan{<:AbstractReconstructionAlgorithm}, name::Symbol, x) = setAll!(plan.parameter, name, x)
function setAll!(plan; kwargs...)
  for key in keys(kwargs)
    setAll!(plan, key, kwargs[key])
  end
end
setAll!(plan::RecoPlan, dict::Dict{Symbol, Any}) = setAll!(plan; dict...)
setAll!(plan::RecoPlan, dict::Dict{String, Any}) = setAll!(plan, Dict{Symbol, Any}(Symbol(k) => v for (k,v) in dict))


export types, type
types(::RecoPlan{T}) where {T<:AbstractReconstructionAlgorithmParameter} = fieldtypes(T)
type(::RecoPlan{T}, name::Symbol) where {T<:AbstractReconstructionAlgorithmParameter} = fieldtype(T, name)

function type(plan::RecoPlan{T}, name::Symbol) where {T<:AbstractReconstructionAlgorithm}
  if name == :parameter
    return RecoPlan
  else
    error("type $(typeof(plan)) has no field $name")
  end
end
types(::RecoPlan{T}) where {T<:AbstractReconstructionAlgorithm} = [type(plan, name) for name in propertynames(plan)]


export clear!
function clear!(plan::RecoPlan{T}, preserve::Bool = true) where {T<:AbstractReconstructionAlgorithmParameter}
  dict = getfield(plan, :values)
  set = getfield(plan, :setProperties)
  for key in keys(dict)
    value = dict[key]
    if typeof(value) <: RecoPlan && preserve 
      clear!(value, preserve)
    else
      dict[key] = missing
      set[key] = false
    end
  end
  return plan
end
clear!(plan::RecoPlan{T}, preserve::Bool = true) where {T<:AbstractReconstructionAlgorithm} = clear!(plan.parameter, preserve)

export getlisteners, addListener!, removeListener!
getlisteners(plan::RecoPlan, field::Symbol) = getfield(plan, :listeners)[field]
function addListener!(plan::RecoPlan, field::Symbol, listener::AbstractPlanListener)
  listeners = getlisteners(plan, field)
  push!(listeners, listener)
end
function removeListener!(plan::RecoPlan, field::Symbol, listener::AbstractPlanListener)
  listeners = getlisteners(plan, field)
  idx = findall(x->isequal(x, listener), listeners)
  isnothing(idx) && deleteat!(listeners, idx)
end

export build
function build(plan::RecoPlan{T}) where {T<:AbstractReconstructionAlgorithmParameter}
  fields = copy(getfield(plan, :values))
  nestedPlans = filter(entry -> isa(last(entry), RecoPlan), fields)
  for (name, nested) in nestedPlans
    fields[name] = build(nested)
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
  plan = RecoPlan(typeof(param))
  for field in fieldnames(typeof(param))
    value = getproperty(param, field)
    if typeof(value) <: AbstractReconstructionAlgorithmParameter || typeof(value) <: AbstractReconstructionAlgorithm
      args[field] = toPlan(plan, value)
    else
      args[field] = value
    end
  end
  setvalues!(plan; args...)
  return plan
end
function toPlan(parent::RecoPlan, x)
  plan = toPlan(x)
  parent!(plan, parent)
  return plan
end 
toPlan(algo::AbstractReconstructionAlgorithm) = toPlan(algo, parameter(algo))
toPlan(algo::AbstractReconstructionAlgorithm, params::AbstractReconstructionAlgorithmParameter) = toPlan(typeof(algo), params) 
function toPlan(::Type{T}, params::AbstractReconstructionAlgorithmParameter) where {T<:AbstractReconstructionAlgorithm}
  plan = RecoPlan(T)
  plan[:parameter] = toPlan(plan, params)
  return plan
end

toDictModule(plan::RecoPlan{T}) where {T} = parentmodule(T)
toDictType(plan::RecoPlan{T}) where {T} = RecoPlan{getfield(parentmodule(T), nameof(T))}
function addDictValue!(dict, value::RecoPlan)
  for field in propertynames(value)
    x = getproperty(value, field)
    if !ismissing(x)
      dict[string(field)] = toDictValue(type(value, field), x)
    end
  end
  listeners = filter(x-> !isempty(last(x)), getfield(value, :listeners))
  if !isempty(listeners)
    listenerDict = Dict{String, Any}()
    for (field, l) in listeners
      serializable = filter(x-> x isa SerializableListener, l)
      if !isempty(serializable)
        listenerDict[string(field)] = toDictValue(typeof(l), l)
      end
    end
    if !isempty(listenerDict) 
      dict[".listener"] = listenerDict
    end
  end
  return dict
end

export plandir, planpath
function plandir(m::Module)
  if m != AbstractImageReconstruction && hasproperty(m, :plandir)
    return getproperty(m, :plandir)()
  else
    return @get_scratch!(string(m))
  end
end
planpath(m::Module, name::AbstractString) = joinpath(plandir(m), string(name, ".toml"))

export loadPlan
loadPlan(m::Module, name::AbstractString, modules::Vector{Module}) = loadPlan(planpath(m, name), modules)
function loadPlan(filename::AbstractString, modules::Vector{Module})
  dict = TOML.parsefile(filename)
  modDict = createModuleDataTypeDict(modules)
  plan = loadPlan!(dict, modDict)
  loadListeners!(plan, dict, modDict)
  return plan
end
function createModuleDataTypeDict(modules::Vector{Module})
  modDict = Dict{String, Dict{String, Union{DataType, UnionAll, Function}}}()
  for mod in modules
    typeDict = Dict{String, Union{DataType, UnionAll, Function}}()
    for field in names(mod)
      try
        t = getfield(mod, field)
        if t isa DataType || t isa UnionAll || t isa Function
          typeDict[string(t)] = t
        end
      catch
      end
    end
    modDict[string(mod)] = typeDict
  end
  return modDict
end
function loadPlan!(dict::Dict{String, Any}, modDict)
  re = r"RecoPlan\{(.*)\}"
  m = match(re, dict[TYPE_TAG])
  if !isnothing(m)
    type = m.captures[1]
    mod = dict[MODULE_TAG]
    plan = RecoPlan(modDict[mod][type])
    loadPlan!(plan, dict, modDict)
    return plan
  else
    # Has to be parameter or algo or broken toml
    # TODO implement
    error("Not implemented yet")
  end
end
function loadPlan!(plan::RecoPlan{T}, dict::Dict{String, Any}, modDict) where {T<:AbstractReconstructionAlgorithm}
  temp = loadPlan!(dict["parameter"], modDict)
  parent!(temp, plan)
  setvalue!(plan, :parameter, temp)
  return plan
end
function loadPlan!(plan::RecoPlan{T}, dict::Dict{String, Any}, modDict) where {T<:AbstractReconstructionAlgorithmParameter}
  for name in propertynames(plan)
    t = type(plan, name)
    param = missing
    key = string(name)
    if haskey(dict, key)
      if t <: AbstractReconstructionAlgorithm || t <: AbstractReconstructionAlgorithmParameter
        param = loadPlan!(dict[key], modDict)
        parent!(param, plan)
      else
        param = loadPlanValue(T, name, t, dict[key], modDict)
      end
    end
    setvalue!(plan, name, param)
  end
  return plan
end
loadPlanValue(parent::Type{T}, field::Symbol, type, value, modDict) where T <: AbstractReconstructionAlgorithmParameter = loadPlanValue(type, value, modDict)
# Type{<:T} where {T}
function loadPlanValue(t::UnionAll, value::Dict, modDict)
  if value[TYPE_TAG] == string(Type)
    return modDict[value[MODULE_TAG]][value[VALUE_TAG]]
  else
    return fromTOML(specializeType(t, value, modDict), value)
  end
end
function loadPlanValue(::Type{Vector{<:T}}, value::Vector, modDict) where {T}
  result = Any[]
  for val in value
    type = modDict[val[MODULE_TAG]][val[TYPE_TAG]]
    push!(result, fromTOML(type, val))
  end
  # Narrow vector
  return identity.(result)
end
uniontypes(t::Union) = Base.uniontypes(t)
#uniontypes(t::Union) = [t.a, uniontypes(t.b)...]
#uniontypes(t::DataType) = [t]
function loadPlanValue(t::Union, value::Dict, modDict)
  types = uniontypes(t)
  idx = findfirst(x-> string(x) == value[UNION_TYPE_TAG], types)
  if isnothing(idx)
    toml = tomlType(value, modDict, prefix = "union")
    idx = !isnothing(toml) ? findfirst(x-> toml <: x, types) : idx # Potentially check if more than one fits and chose most fitting
  end
  type = isnothing(idx) ? t : types[idx]
  return loadPlanValue(type, value[VALUE_TAG], modDict)
end
loadPlanValue(t::DataType, value::Dict, modDict) = fromTOML(specializeType(t, value, modDict), value)
loadPlanValue(t, value, modDict) = fromTOML(t, value)

function tomlType(dict::Dict, modDict; prefix::String = "")
  if haskey(dict, ".$(prefix)module") && haskey(dict, ".$(prefix)type")
    mod = dict[".$(prefix)module"]
    type = dict[".$(prefix)type"]
    if haskey(modDict, mod) && haskey(modDict[mod], type)
      return modDict[mod][type]
    end
  end
  return nothing
end
function specializeType(t::Union{DataType, UnionAll}, value::Dict, modDict)
  if isconcretetype(t)
    return t
  end
  type = tomlType(value, modDict)
  return !isnothing(type) && type <: t ? type : t 
end

loadListeners!(plan, dict, modDict) = loadListeners!(plan, plan, dict, modDict)
function loadListeners!(root::RecoPlan, plan::RecoPlan{T}, dict, modDict) where {T<:AbstractReconstructionAlgorithm}
  loadListeners!(root, plan.parameter, dict["parameter"], modDict)
end
function loadListeners!(root::RecoPlan, plan::RecoPlan{T}, dict, modDict) where {T<:AbstractReconstructionAlgorithmParameter}
  if haskey(dict, ".listener")
    for (property, listenerDicts) in dict[".listener"]
      for listenerDict in listenerDicts
        listener = loadListener(root, listenerDict, modDict)
        addListener!(plan, Symbol(property), listener)
      end
    end
  end
  for property in propertynames(plan)
    value = plan[property]
    if value isa RecoPlan
      loadListeners!(root, value, dict[string(property)], modDict)
    end
  end
end
export loadListener
function loadListener(root, dict, modDict)
  type = tomlType(dict, modDict)
  return loadListener(type, root, dict, modDict)
end

export savePlan
savePlan(filename::AbstractString, plan::RecoPlan) = toTOML(filename, plan)
savePlan(m::Module, planname::AbstractString, plan::RecoPlan) = savePlan(planpath(m, planname), plan)

include("LinkedFieldListener.jl")