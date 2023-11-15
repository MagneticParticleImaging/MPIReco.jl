export SystemMatrixWrapper, rows, rows!, parentrows, frequencies!

mutable struct SystemMatrixWrapper{T, P<:AbstractArray{T, 2}} <: AbstractArray{T, 2}
  parent::P
  parentrows::Vector{Int64}
  rows::Vector{Int64}
  smView::SubArray{T, 2, P}
end

# TODO Get original rows for SM somehow
SystemMatrixWrapper(sm::AbstractArray, rows::Vector{Int64}) = SystemMatrixWrapper(sm, rows, rows, view(sm, :, :))
SystemMatrixWrapper(sm::AbstractArray, freqs::Vector{Int64}, freqsize::Int64) = SystemMatrixWrapper(sm, frequenciesToRows(freqs, freqsize))


Base.parent(sm::SystemMatrixWrapper) = sm.parent
parentrows(sm::SystemMatrixWrapper) = sm.parentrows

# maybe not rows! but components?
function rows!(sm::SystemMatrixWrapper, rows::Vector{Int64})
  if !issubset(rows, parentrows(sm))
    error("Given rows are not a subset of parent rows")
  end
  mapping = map(x-> findfirst(y-> x==y, sm.parentrows), rows)
  if sm.parent isa Transpose
    sm.smView = view(sm.parent, mapping, :)
  else
    sm.smView = view(sm.parent, :, mapping)
  end
  return sm
end
rows(sm::SystemMatrixWrapper) = sm.rows
frequencies!(sm::SystemMatrixWrapper, freqs::Vector{Int64}, freqsize::Int64) = rows!(sm, frequenciesToRows(freqs, freqsize))

frequenciesToRows(freqs::Vector{Int64}, freqsize::Int64) = vcat(collect.(map(x-> (x-1)*freqsize+1:x*freqsize, freqs))...)

Base.getindex(sm::SystemMatrixWrapper, i::Vararg{Int, N}) where {N} = sm.smView[i...]
Base.setindex!(sm::SystemMatrixWrapper, x, i::Vararg{Int, N}) where {N} = sm.smView[i...] = x
Base.firstindex(sm::SystemMatrixWrapper) = firstindex(sm.smView)
Base.lastindex(sm::SystemMatrixWrapper) = lastindex(sm.smView)
Base.size(sm::SystemMatrixWrapper) = size(sm.smView)
Base.length(sm::SystemMatrixWrapper) = length(sm.smView)

# TODO
#dot_with_matrix_row
#mul!
#mul! with adjoint(SystemMatrixWrapper)