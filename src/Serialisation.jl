RecoUtils.toDictValue(file::MPIFile) = filepath(file)
RecoUtils.toDictValue(file::Union{MultiMPIFile, MultiContrastFile}) = filepath.(file)
RecoUtils.toDictValue(norm::AbstractRegularizationNormalization) = RecoUtils.toDict(norm)

# TODO Check how nothing is stored and maybe return it here
function RecoUtils.fromTOML(::Type{T}, dict::Dict{String, Any}) where {T<: UnitRange}
  start = dict["start"]
  stop = dict["stop"]
  return UnitRange(start, stop)
end
RecoUtils.fromTOML(::Type{<:MPIFile}, x::AbstractString) = MPIFile(x)
RecoUtils.fromTOML(::Type{<:MultiMPIFile}, x::Vector) = MultiMPIFile(MPIFile.(x))
function RecoUtils.fromTOML(::Type{T}, dict::Dict{String, Any}) where {T<: AbstractRegularization}
  λ = dict["λ"]
  filteredKeys = filter(!isequal("λ"), keys(dict))
  kwargs = Dict(Symbol(x) => dict[x] for x in filteredKeys)
  return T(λ; kwargs...)  
end