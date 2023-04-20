export toTOML, toDict, toDict!, toDictValue

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

function toKwargs(value; flatten::Vector{DataType} = DataType[])
  dict = Dict{Symbol, Any}()
  return toKwargs!(dict, value, flatten = flatten)
end

function toKwargs(values::Vector; flatten::Vector{DataType} = DataType[])
  dict = Dict{Symbol, Any}()
  foreach(i-> toKwargs!(dict, i, flatten = flatten), values)
  return dict
end

function toKwargs!(dict, value; flatten::Vector{DataType} = DataType[])
  for field in fieldnames(typeof(value))
    prop = getproperty(value, field)
    if any(i -> prop isa i, flatten)
      toKwargs!(dict, prop, flatten = flatten)
    else
      dict[field] = prop
    end
  end
  return dict
end

function toKwargs(v::Union{AbstractReconstructionAlgorithmParameter, Vector{AbstractReconstructionAlgorithmParameter}})
  return toKwargs(v, flatten = [AbstractReconstructionAlgorithmParameter])
end

function fromKwargs(type::Type{T}, kargs...) where {T}
  args = Dict{Symbol, Any}()
  dict = values(kargs)
  for field in fieldnames(type)
    if haskey(dict, field)
      args[field] = dict.field
    end
  end
  return type(args)
end