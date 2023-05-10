export toTOML, toDict, toDict!, toDictValue, toKwargs, toKwargs!, fromKwargs

function toTOML(fileName::AbstractString, value)
  open(fileName, "w") do io
    toTOML(io, value)
  end
end

function toTOML(io::IO, value)
  dict = toDict(value)
  TOML.print(io, dict)
end

function toDict(value)
  dict = Dict{String, Any}()
  return toDict!(dict, value)
end

function toDict!(dict, value)
  for field in fieldnames(typeof(value))
    dict[String(field)] = toDictValue(getproperty(value, field))
  end
  return dict
end

function toDictValue(x)
  if fieldcount(typeof(x)) > 0
    return toDict(x)
  else
    return x
  end
end
toDictValue(x::T) where {T<:Enum} = string(x)
toDictValue(x::Array) = toDictValue.(x)

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