function Adapt.adapt_structure(::Type{arrT}, op::MultiPatchOperator) where {arrT <: AbstractGPUArray}
  validSMs = all(x -> size(x) == size(op.S[1]), op.S)
  validXCC = all(x -> length(x) == length(op.xcc[1]), op.xcc)
  validXSS = all(x -> length(x) == length(op.xss[1]), op.xss)

  # Ideally we create a DenseMultiPatchOperator on the GPU
  if validSMs && validXCC && validXSS
    S = adapt(arrT, stack(op.S))
    # We want to use Int32 for better GPU performance
    xcc = Int32.(adapt(arrT, stack(op.xcc)))
    xss = Int32.(adapt(arrT, stack(op.xss)))
    sign = Int32.(adapt(arrT, op.sign))
    RowToPatch = Int32.(adapt(arrT, op.RowToPatch))
    patchToSMIdx = Int32.(adapt(arrT, op.patchToSMIdx))
    return DenseMultiPatchOperator(S, op.grid, Int32(op.N), Int32(op.M), RowToPatch, xcc, xss, sign, Int32(op.nPatches), patchToSMIdx)
  else
    throw(ArgumentError("Cannot adapt MultiPatchOperator to $arrT, since it cannot be represented as a DenseMultiPatchOperator"))
  end
end

@kernel function dense_mul!(b, @Const(x), @Const(S), @Const(xcc), @Const(xss), @Const(signs), @Const(M), @Const(RowToPatch), @Const(patchToSMIdx))
  # Each group/block handles a single row of the operator
  operator_row = @index(Group, Linear) # k
  patch = RowToPatch[operator_row] # p
  patch_row = mod1(operator_row, M) # j
  smIdx = patchToSMIdx[patch]
  sign = signs[patch_row, smIdx]
  grid_stride = prod(@groupsize())
  N = size(xss, 1)
  
  # We want to use a grid-stride loop to perform the sparse matrix-vector product.
  # Each thread performs a single element-wise multiplication and reduction in its shared spot.
  # Afterwards we reduce over the shared memory.
  localIdx = @index(Local, Linear)
  shared = @localmem eltype(b) prod(@groupsize())
  shared[localIdx] = zero(eltype(b))

  # First we iterate over the sparse indices
  i = localIdx
  while i <= N
    shared[localIdx] = shared[localIdx] + sign * S[patch_row, xss[i, patch], smIdx] * x[xcc[i, patch]]
    i += grid_stride
  end
  @synchronize

  # Now we need to reduce the shared memory to get the final result
  @private s = div(min(grid_stride, N), Int32(2))
  while s > Int32(0)
    if localIdx <= s
      shared[localIdx] = shared[localIdx] + shared[localIdx + s]
    end
    s >>= 1
    @synchronize
  end

  # Write the result out to b
  if localIdx == 1
    b[operator_row] = shared[localIdx]
  end
end

function LinearAlgebra.mul!(b::AbstractVector{T}, op::DenseMultiPatchOperator{T, V}, x::AbstractVector{T}) where {T, V}
  b[:] .= zero(T)
  backend = get_backend(b)
  kernel = dense_mul!(backend, 256)
  kernel(b, x, op.S, op.xcc, op.xss, op.sign, div(op.M, op.nPatches), op.RowToPatch, op.patchToSMIdx; ndrange = (256, size(op, 1)))
  synchronize(backend)
  return b
end