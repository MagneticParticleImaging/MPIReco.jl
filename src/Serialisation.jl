AbstractImageReconstruction.toDictValue(file::MPIFile) = AbstractImageReconstruction.toDict(file)
AbstractImageReconstruction.toDictValue!(dict, file::Union{MultiMPIFile, MultiContrastFile}) = dict[AbstractImageReconstruction.VALUE_TAG] = filepath.(file)
AbstractImageReconstruction.toDictValue!(dict, file::MPIFile) = dict[AbstractImageReconstruction.VALUE_TAG] = filepath(file)
AbstractImageReconstruction.toDictValue(norm::AbstractRegularizationNormalization) = AbstractImageReconstruction.toDict(norm)

function AbstractImageReconstruction.fromTOML(::Type{T}, dict::Dict{String, Any}) where {T<: UnitRange}
  start = dict["start"]
  stop = dict["stop"]
  return UnitRange(start, stop)
end
AbstractImageReconstruction.fromTOML(t::Type{<:MPIFile}, x::Dict) = AbstractImageReconstruction.fromTOML(t, x[AbstractImageReconstruction.VALUE_TAG])
AbstractImageReconstruction.fromTOML(::Type{<:MPIFile}, x::AbstractString) = MPIFile(x)
AbstractImageReconstruction.fromTOML(::Type{<:MultiMPIFile}, x::Vector) = MultiMPIFile(MPIFile.(x))
function AbstractImageReconstruction.fromTOML(::Type{T}, dict::Dict{String, Any}) where {T<: AbstractRegularization}
  λ = dict["λ"]
  filteredKeys = filter(x-> !(isequal("λ",x) || startswith(x, ".")), keys(dict))
  kwargs = Dict(Symbol(x) => dict[x] for x in filteredKeys)
  if haskey(kwargs, :sparseTrafo) && isempty(kwargs[:sparseTrafo])
    pop!(kwargs, :sparseTrafo)
  end
  return T(λ; kwargs...)  
end
AbstractImageReconstruction.fromTOML(t::Type{T}, dict) where {T<:AbstractRegularizationNormalization} = T()