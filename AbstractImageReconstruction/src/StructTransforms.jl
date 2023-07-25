export toTOML, toDict, toDict!, toDictValue, toKwargs, toKwargs!, fromKwargs

# TODO adapt tomlType
const MODULE_TAG = ".module"
const TYPE_TAG = ".type"
const VALUE_TAG = ".value"
const UNION_TYPE_TAG = ".uniontype"
const UNION_MODULE_TAG = ".unionmodule"

function toTOML(fileName::AbstractString, value)
  open(fileName, "w") do io
    toTOML(io, value)
  end
end

function toTOML(io::IO, value)
  dict = toDict(value)
  TOML.print(io, dict) do x
    toTOML(x)
  end
end

toTOML(x::Module) = string(x)
toTOML(x::Symbol) = string(x)
toTOML(x::T) where {T<:Enum} = string(x)
toTOML(x::Array) = toTOML.(x)
toTOML(x::Type{T}) where T = string(x)
toTOML(x::Nothing) = Dict()

fromTOML(t, x) = x
function fromTOML(::Type{Nothing}, x::Dict) where {T}
  if isempty(x)
    return nothing
  end
  error("Unexpected value $x for Nothing, expected empty Dict")
end
fromTOML(::Type{V}, x::Vector) where {T, V<:Vector{<:T}} = fromTOML.(T, x)

function toDict(value)
  dict = Dict{String, Any}()
  return toDict!(dict, value)
end

function toDict!(dict, value)
  dict[MODULE_TAG] = toDictModule(value)
  dict[TYPE_TAG] = toDictType(value)
  addDictValue!(dict, value)
  return dict
end
toDictModule(value) = parentmodule(typeof(value))
toDictType(value) = nameof(typeof(value))
function addDictValue!(dict, value)
  for field in fieldnames(typeof(value))
    dict[string(field)] = toDictValue(fieldtype(typeof(value), field), getfield(value, field))
  end
end

toDictValue(type, value) = toDictValue(value)
function toDictValue(x)
  if fieldcount(typeof(x)) > 0
    return toDict(x)
  else
    return x
  end
end
toDictValue(x::Array) = toDictValue.(x)
toDictValue(x::Type{T}) where T = toDict(x)
function toDict!(dict, ::Type{T}) where T
  dict[MODULE_TAG] = parentmodule(T)
  dict[TYPE_TAG] = Type
  dict[VALUE_TAG] = T
  return dict
end

function toDictValue(type::Union, value)
  dict = Dict{String, Any}()
  dict[MODULE_TAG] = toDictModule(type)
  dict[TYPE_TAG] = toDictType(type)
  dict[VALUE_TAG] = toDictValue(value)
  dict[UNION_TYPE_TAG] = typeof(value) # directly type to not remove parametric fields
  dict[UNION_MODULE_TAG] = toDictModule(typeof(value))
  return dict
end

function toKwargs(value; kwargs...)
  dict = Dict{Symbol, Any}()
  return toKwargs!(dict, value; kwargs...)
end

function toKwargs(values::Vector; kwargs...)
  dict = Dict{Symbol, Any}()
  foreach(i-> toKwargs!(dict, i; kwargs...), values)
  return dict
end

function toKwargs!(dict, value; flatten::Vector{DataType} = DataType[], ignore::Vector{Symbol} = Symbol[], default::Dict{Symbol, Any} = Dict{Symbol, Any}(), overwrite::Dict{Symbol, Any} = Dict{Symbol, Any}())
  for field in fieldnames(typeof(value))
    prop = getproperty(value, field)
    if in(field, ignore)
      # NOP
    elseif any(i -> prop isa i, flatten)
      toKwargs!(dict, prop, flatten = flatten, ignore = ignore, default = default)
    elseif (isnothing(prop) || ismissing(prop)) && haskey(default, field)
      dict[field] = default[field]
    else
      dict[field] = prop
    end
  end
  for key in keys(overwrite)
    dict[key] = overwrite[key] 
  end
  return dict
end

function toKwargs(v::AbstractReconstructionAlgorithmParameter; flatten::Union{Vector{DataType}, Nothing} = nothing, kwargs...)
  dict = Dict{Symbol, Any}()
  return toKwargs!(dict, v; flatten = isnothing(flatten) ? [AbstractReconstructionAlgorithmParameter] : flatten, kwargs...)
end
function toKwargs(v::Vector{<:AbstractReconstructionAlgorithmParameter}; flatten::Union{Vector{DataType}, Nothing} = nothing, kwargs...)
  dict = Dict{Symbol, Any}()
  flatten = isnothing(flatten) ? [AbstractReconstructionAlgorithmParameter] : flatten
  foreach(i-> toKwargs!(dict, i; flatten = flatten, kwargs...), v)
  return dict
end


function fromKwargs(type::Type{T}; kwargs...) where {T}
  args = Dict{Symbol, Any}()
  dict = values(kwargs)
  for field in fieldnames(type)
    if haskey(dict, field)
      args[field] = getproperty(dict, field)
    end
  end
  return type(;args...)
end