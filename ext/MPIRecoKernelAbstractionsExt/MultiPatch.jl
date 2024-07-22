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
    return DenseMultiPatchOperator(S, op.grid, op.N, op.M, RowToPatch, xcc, xss, sign, Int32(op.nPatches), patchToSMIdx)
  else
    throw(ArgumentError("Cannot adapt MultiPatchOperator to $arrT, since it cannot be represented as a DenseMultiPatchOperator"))
  end
end

@kernel cpu = false inbounds = true function dense_mul!(b, @Const(x), @Const(S), @Const(xcc), @Const(xss), @Const(signs), @Const(M), @Const(RowToPatch), @Const(patchToSMIdx))
  # Each group/block handles a single row of the operator
  operator_row = @index(Group, Linear) # k
  patch = RowToPatch[operator_row] # p
  patch_row = mod1(operator_row, M) # j
  smIdx = patchToSMIdx[patch]
  sign = eltype(b)(signs[patch_row, smIdx])
  grid_stride = prod(@groupsize())
  N = Int32(size(xss, 1))
  
  # We want to use a grid-stride loop to perform the sparse matrix-vector product.
  # Each thread performs a single element-wise multiplication and reduction in its shared spot.
  # Afterwards we reduce over the shared memory.
  localIdx = @index(Local, Linear)
  shared = @localmem eltype(b) grid_stride
  shared[localIdx] = zero(eltype(b))

  # First we iterate over the sparse indices
  @unroll for i = localIdx:grid_stride:N
    shared[localIdx] = shared[localIdx] + sign * S[patch_row, xss[i, patch], smIdx] * x[xcc[i, patch]]
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

function LinearAlgebra.mul!(b::AbstractVector{T}, op::DenseMultiPatchOperator{T, V}, x::AbstractVector{T}) where {T, V <: AbstractGPUArray}
  backend = get_backend(b)
  kernel = dense_mul!(backend, 256)
  kernel(b, x, op.S, op.xcc, op.xss, op.sign, Int32(div(op.M, op.nPatches)), op.RowToPatch, op.patchToSMIdx; ndrange = (256, size(op, 1)))
  synchronize(backend)
  return b
end

@kernel inbounds = true function dense_mul_adj!(res, @Const(t), @Const(S), @Const(xcc), @Const(xss), @Const(signs), @Const(M), @Const(RowToPatch), @Const(patchToSMIdx))
  # Each group/block handles a single column of the adjoint(operator)
  # i.e. a row of the operator
  localIdx = @index(Local, Linear)
  groupIdx = @index(Group, Linear) # k
  patch = RowToPatch[groupIdx] # p
  patch_row = mod1(groupIdx, M) # j
  smIdx = patchToSMIdx[patch]
  sign = eltype(res)(signs[patch_row, smIdx])
  grid_stride = prod(@groupsize())
  N = Int32(size(xss, 1))
  
  
  # Each thread within the block will add the same value of t
  val = t[groupIdx]
 
  # Since we go along the columns during a matrix-vector product,
  # we have a race condition with other threads writing to the same result.
  for i = localIdx:grid_stride:N
    tmp = sign * adjoint(S[patch_row, xss[i, patch], smIdx]) * val
    # @atomic is not supported for ComplexF32 numbers
    Atomix.@atomic res[1, xcc[i, patch]] += tmp.re
    Atomix.@atomic res[2, xcc[i, patch]] += tmp.im
  end
end

function LinearAlgebra.mul!(res::AbstractVector{T}, adj::Adjoint{T, OP}, t::AbstractVector{T}) where {T <: Complex, V <: AbstractArray, OP <: DenseMultiPatchOperator{T, V}}
  backend = get_backend(res)
  op = adj.parent
  res .= zero(T) # We need to zero the result, because we are using += in the kernel
  kernel = dense_mul_adj!(backend, 256)
  # We have to reinterpret the result as a real array, because atomic operations on Complex numbers are not supported
  kernel(reinterpret(reshape, real(eltype(res)), res), t, op.S, op.xcc, op.xss, op.sign, Int32(div(op.M, op.nPatches)), op.RowToPatch, op.patchToSMIdx; ndrange = (256, size(op, 1)))
  synchronize(backend)
  return res
end