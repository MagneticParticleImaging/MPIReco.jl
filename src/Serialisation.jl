RecoUtils.toDictValue(file::MPIFile) = RecoUtils.toDict(file)
RecoUtils.addDictValue!(dict, file::Union{MultiMPIFile, MultiContrastFile}) = dict[RecoUtils.VALUE_TAG] = filepath.(file)
RecoUtils.addDictValue!(dict, file::MPIFile) = dict[RecoUtils.VALUE_TAG] = filepath(file)
RecoUtils.toDictValue(norm::AbstractRegularizationNormalization) = RecoUtils.toDict(norm)

function RecoUtils.fromTOML(::Type{T}, dict::Dict{String, Any}) where {T<: UnitRange}
  start = dict["start"]
  stop = dict["stop"]
  return UnitRange(start, stop)
end
RecoUtils.fromTOML(t::Type{<:MPIFile}, x::Dict) = RecoUtils.fromTOML(t, x[RecoUtils.VALUE_TAG])
RecoUtils.fromTOML(::Type{<:MPIFile}, x::AbstractString) = MPIFile(x)
RecoUtils.fromTOML(::Type{<:MultiMPIFile}, x::Vector) = MultiMPIFile(MPIFile.(x))
function RecoUtils.fromTOML(::Type{T}, dict::Dict{String, Any}) where {T<: AbstractRegularization}
  λ = dict["λ"]
  filteredKeys = filter(x-> !(isequal("λ",x) || startswith(x, ".")), keys(dict))
  kwargs = Dict(Symbol(x) => dict[x] for x in filteredKeys)
  if haskey(kwargs, :sparseTrafo) && isempty(kwargs[:sparseTrafo])
    pop!(kwargs, :sparseTrafo)
  end
  return T(λ; kwargs...)  
end
RecoUtils.fromTOML(t::Type{T}, dict) where {T<:AbstractRegularizationNormalization} = T()