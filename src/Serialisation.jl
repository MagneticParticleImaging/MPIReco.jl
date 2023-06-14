RecoUtils.toDictValue(file::MPIFile) = filepath(file)
RecoUtils.toDictValue(norm::AbstractRegularizationNormalization) = RecoUtils.toDict(norm)

# TODO Check how nothing is stored and maybe return it here
RecoUtils.fromTOML(::Type{Union{Nothing, T}}, x) where {T} = RecoUtils.fromTOML(T, x)
function RecoUtils.fromTOML(::Type{<:UnitRange}, dict::Dict{String, Any})
  start = dict["start"]
  stop = dict["stop"]
  return UnitRange(start, stop)
end
RecoUtils.fromTOML(::Type{<:MPIFile}, x::AbstractString) = MPIFile(x)
function RecoUtils.fromTOML(::Type{T}, dict::Dict{String, Any}) where {T<: AbstractRegularization}
  λ = dict["λ"]
  filteredKeys = filter(!isequal("λ"), keys(dict))
  kwargs = Dict(Symbol(x) => dict[x] for x in filteredKeys)
  return T(λ; kwargs...)  
end