function Adapt.adapt_structure(::Type{arrT}, op::MultiPatchOperator) where {arrT <: AbstractGPUArray}
  validSMs = all(x -> size(x) == size(op.S[1]), op.S)
  validXCC = all(x -> length(x) == length(op.xcc[1]), op.xcc)
  validXSS = all(x -> length(x) == length(op.xss[1]), op.xss)

  # Ideally we create a DenseMultiPatchOperator on the GPU
  if validSMs && validXCC && validXSS
    # We want to use Int32 for better GPU performance
    xcc = Int32.(stack(adapt.(arrT, op.xcc)))
    xss = Int32.(stack(adapt.(arrT, op.xss)))
    sign = Int32.(adapt(arrT, op.sign))
    RowToPatch = Int32.(adapt(arrT, op.RowToPatch))
    patchToSMIdx = Int32.(adapt(arrT, op.patchToSMIdx))

    backend = KernelAbstractions.get_backend(xcc)
    S = KernelAbstractions.allocate(backend, eltype(first(op.S)), size(transpose(first(op.S)))..., length(op.S))
    cache = GPUArrays.AllocCache()
    for (i, sm) in enumerate(op.S)
      GPUArrays.@cached cache S[:, :, i] .= transpose(adapt(arrT, sm))
    end
    GPUArrays.unsafe_free!(cache)

    return DenseMultiPatchOperator(S, op.grid, op.N, op.M, RowToPatch, xcc, xss, sign, Int32(op.nPatches), patchToSMIdx)
  else
    throw(ArgumentError("Cannot adapt MultiPatchOperator to $arrT, since it cannot be represented as a DenseMultiPatchOperator"))
  end
end


function LinearAlgebra.mul!(b::AbstractVector{T}, op::DenseMultiPatchOperator{T, V}, x::AbstractVector{T}) where {T, V <: AbstractGPUArray}
  backend = get_backend(b)


  # Define kernel within the mul! function to make sure we know the workgroup size for manual unrolling
  @kernel cpu = false inbounds = true function dense_mul!(b, @Const(x), @Const(S), @Const(xcc), @Const(xss), @Const(signs), @Const(M), @Const(RowToPatch), @Const(patchToSMIdx))
    # Each group/block handles a single row of the operator
    operator_row = @index(Group, Linear) # k
    patch = RowToPatch[operator_row] # p
    patch_row = mod1(operator_row, M) # j
    smIdx = patchToSMIdx[patch]
    sign = eltype(b)(signs[patch_row, smIdx])
    @uniform grid_stride = prod(@groupsize())
    N = Int32(size(xss, 1))
    
    # We want to use a grid-stride loop to perform the sparse matrix-vector product.
    # Each thread performs a single element-wise multiplication and reduction in its shared spot.
    # Afterwards we reduce over the shared memory.
    localIdx = @index(Local, Linear)
    shared = @localmem eltype(b) grid_stride
    shared[localIdx] = zero(eltype(b))
  
    # First we iterate over the sparse indices
    tmp = zero(eltype(b))
    @unroll for i = localIdx:grid_stride:N
      tmp += sign * S[xss[i, patch], patch_row, smIdx] * x[xcc[i, patch]]
    end
    shared[localIdx] = tmp
    # We first sum in a temp variable, hoping that it is accumulated in a register, since registers are faster than shared memory
    @synchronize
  
    # Now we need to reduce the shared memory to get the final result
    full_reduction = grid_stride < N
    if full_reduction
      
      # For a full reduction we know s = 512 and can (manually) unroll our loop
      #localIdx <= 512 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx + 512])
      #@synchronize
      localIdx <= 256 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx + 256])
      @synchronize
      localIdx <= 128 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx + 128])
      @synchronize
      localIdx <= 64 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx + 64])
      @synchronize
      localIdx <= 32 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx + 32])
      @synchronize
      localIdx <= 16 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx + 16])
      @synchronize
      localIdx <= 8 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx + 8])
      @synchronize
      localIdx <= 4 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx + 4])
      @synchronize
      localIdx <= 2 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx + 2])
      @synchronize
      localIdx == 1 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx + 1])
      
   
    else
      @private s = div(min(grid_stride, N), Int32(2))
      while s > Int32(0)
        if localIdx <= s
          shared[localIdx] = shared[localIdx] + shared[localIdx + s]
        end
        s >>= 1
        @synchronize
      end
    end
  
    # Write the result out to b
    if localIdx == 1
      b[operator_row] = shared[localIdx]
    end
  end

  kernel = dense_mul!(backend, 512, (512, size(op, 1)))
  kernel(b, x, op.S, op.xcc, op.xss, op.sign, Int32(div(op.M, op.nPatches)), op.RowToPatch, op.patchToSMIdx; ndrange = (512, size(op, 1)))
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
  @uniform grid_stride = prod(@groupsize())
  N = Int32(size(xss, 1))
  
  
  # Each thread within the block will add the same value of t
  val = t[groupIdx]
 
  # Since we go along the columns during a matrix-vector product,
  # we have a race condition with other threads writing to the same result.
  @unroll for i = localIdx:grid_stride:N
    tmp = sign * conj(S[xss[i, patch], patch_row, smIdx]) * val
    # @atomic is not supported for ComplexF32 numbers
    Atomix.@atomic res[1, xcc[i, patch]] += real(tmp)
    Atomix.@atomic res[2, xcc[i, patch]] += imag(tmp)
  end
end

function LinearAlgebra.mul!(res::AbstractVector{T}, adj::Adjoint{T, OP}, t::AbstractVector{T}) where {T <: Complex, V <: AbstractGPUArray, OP <: DenseMultiPatchOperator{T, V}}
  backend = get_backend(res)
  op = adj.parent
  res .= zero(T) # We need to zero the result, because we are using += in the kernel
  kernel = dense_mul_adj!(backend, 512, (512, size(op, 1)))
  # We have to reinterpret the result as a real array, because atomic operations on Complex numbers are not supported
  kernel(reinterpret(reshape, real(eltype(res)), res), t, op.S, op.xcc, op.xss, op.sign, Int32(div(op.M, op.nPatches)), op.RowToPatch, op.patchToSMIdx; ndrange = (512, size(op, 1)))
  return res
end

function LinearAlgebra.mul!(res::AbstractVector{Tc}, nop::NormalOp{Tc, OP}, x::AbstractVector) where {T, Tc <: Complex{T}, V <: AbstractGPUArray, OP <: DenseMultiPatchOperator{Tc, V}}
  weights = prepareKernelWeights(T, nop.weights)
  return mul_dense_normal!(res, nop, x, weights)
end
# Known weights
prepareKernelWeights(T, weights::WeightingOp) = weights.weights
prepareKernelWeights(::Type{T}, weights::Nothing) where T= one(T)
# Unknown weight, cant do kernel fusion
prepareKernelWeights(T, weights) = nothing

function mul_dense_normal!(res::AbstractVector{Tc}, nop::NormalOp{Tc, OP}, x::AbstractVector, weights::Nothing) where {T, Tc <: Complex{T}, V <: AbstractGPUArray, OP <: DenseMultiPatchOperator{Tc, V}}
  op = nop.parent
  mul!(nop.tmp, op, x)
  mul!(nop.tmp, nop.weights, nop.tmp)
  mul!(res, adjoint(op), nop.tmp)
end

function mul_dense_normal!(res::AbstractVector{Tc}, nop::NormalOp{Tc, OP}, x::AbstractVector, weights) where {T, Tc <: Complex{T}, V <: AbstractGPUArray, OP <: DenseMultiPatchOperator{Tc, V}}
  backend = get_backend(res)
  op = nop.parent
  res .= zero(T) # We need to zero the result, because we are using += in the kernel

  @kernel cpu = false function dense_mul_normal!(res, @Const(x), @Const(S), @Const(xcc), @Const(xss), @Const(signs), @Const(M), @Const(RowToPatch), @Const(patchToSMIdx), @Const(weights))
    ### Forward operator ###
    # Each group/block handles a single row of the operator
    operator_row = @index(Group, Linear) # k
    patch = RowToPatch[operator_row] # p
    patch_row = mod1(operator_row, M) # j
    smIdx = patchToSMIdx[patch]
    sign = eltype(x)(signs[patch_row, smIdx])
    @uniform grid_stride = prod(@groupsize())
    N = Int32(size(xss, 1))

    # We want to use a grid-stride loop to perform the sparse matrix-vector product.
    # Each thread performs a single element-wise multiplication and reduction in its shared spot.
    # Afterwards we reduce over the shared memory.
    localIdx = @index(Local, Linear)
    shared = @localmem eltype(x) grid_stride
    shared[localIdx] = zero(eltype(x))

    # First we iterate over the sparse indices
    tmp = zero(eltype(x))
    @unroll for i = localIdx:grid_stride:N
      tmp += sign * S[xss[i, patch], patch_row, smIdx] * x[xcc[i, patch]]
    end
    shared[localIdx] = tmp
    # We first sum in a temp variable, hoping that it is accumulated in a register, since registers are faster than shared memory
    @synchronize

    # Now we need to reduce the shared memory to get the final result
    full_reduction = grid_stride < N
    if full_reduction

      # For a full reduction we know s = 512 and can (manually) unroll our loop
      #localIdx <= 512 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx + 512])
      #@synchronize
      localIdx <= 256 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx+256])
      @synchronize
      localIdx <= 128 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx+128])
      @synchronize
      localIdx <= 64 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx+64])
      @synchronize
      localIdx <= 32 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx+32])
      @synchronize
      localIdx <= 16 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx+16])
      @synchronize
      localIdx <= 8 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx+8])
      @synchronize
      localIdx <= 4 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx+4])
      @synchronize
      localIdx <= 2 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx+2])
      @synchronize
      localIdx == 1 && (@inbounds shared[localIdx] = shared[localIdx] + shared[localIdx+1])
      @synchronize


    else
      @private s = div(min(grid_stride, N), Int32(2))
      while s > Int32(0)
        if localIdx <= s
          shared[localIdx] = shared[localIdx] + shared[localIdx+s]
        end
        s >>= 1
        @synchronize
      end
    end

    ### Adjoint operator ###
    val = shared[1] * get_kernel_weights(weights, operator_row)
    @unroll for i = localIdx:grid_stride:N
      tmp2 = sign * conj(S[xss[i, patch], patch_row, smIdx]) * val
      # @atomic is not supported for ComplexF32 numbers
      Atomix.@atomic res[1, xcc[i, patch]] += real(tmp2)
      Atomix.@atomic res[2, xcc[i, patch]] += imag(tmp2)
    end  
  end

  kernel = dense_mul_normal!(backend, 512, (512, size(op, 1)))
  kernel(reinterpret(reshape, T, res), x, op.S, op.xcc, op.xss, op.sign, Int32(div(op.M, op.nPatches)), op.RowToPatch, op.patchToSMIdx, weights; ndrange = (512, size(op, 1)))
  return res
end

# Kaczmarz specific functions
function RegularizedLeastSquares.dot_with_matrix_row(op::DenseMultiPatchOperator{T, V}, x::AbstractArray{T}, k::Int) where {T, V <: AbstractGPUArray}
  patch = @allowscalar op.RowToPatch[k]
  patch_row = mod1(k, div(op.M,op.nPatches))
  smIdx = @allowscalar op.patchToSMIdx[patch]
  sign = @allowscalar op.sign[patch_row, smIdx]
  S = op.S
  # Inplace reduce-broadcast: https://github.com/JuliaLang/julia/pull/31020
  return sum(Broadcast.instantiate(Base.broadcasted(view(op.xss, :, patch), view(op.xcc, :, patch)) do xs, xc
    @inbounds sign * S[xs, patch_row, smIdx] * x[xc]
  end))
end

function RegularizedLeastSquares.rownorm²(op::DenseMultiPatchOperator{T, V}, row::Int64) where {T, V <: AbstractGPUArray}
  patch = @allowscalar op.RowToPatch[row]
  patch_row = mod1(row, div(op.M,op.nPatches))
  smIdx = @allowscalar op.patchToSMIdx[patch]
  sign = @allowscalar op.sign[patch_row, smIdx]
  S = op.S
  return mapreduce(xs -> abs2(sign * S[xs, patch_row, smIdx]), +, view(op.xss, :, patch))
end

@kernel cpu = false function kaczmarz_update_kernel!(x, @Const(S), @Const(row), @Const(beta), @Const(xcc), @Const(xss), @Const(signs), @Const(M), @Const(RowToPatch), @Const(patchToSMIdx))
  # Each thread handles one element of the kaczmarz update
  idx = @index(Global, Linear)
  patch = RowToPatch[row]
  patch_row = mod1(row, M)
  smIdx = patchToSMIdx[patch]
  sign = eltype(x)(signs[patch_row, smIdx])
  x[xcc[idx, patch]] += beta * conj(sign * S[xss[idx, patch], patch_row, smIdx])
end

function RegularizedLeastSquares.kaczmarz_update!(op::DenseMultiPatchOperator{T, V}, x::vecT, row, beta) where {T, vecT <: AbstractGPUVector{T}, V <: AbstractGPUArray{T}}
  backend = get_backend(x)
  kernel = kaczmarz_update_kernel!(backend, 512)
  kernel(x, op.S, row, beta, op.xcc, op.xss, op.sign, Int32(div(op.M, op.nPatches)), op.RowToPatch, op.patchToSMIdx; ndrange = size(op.xss, 1))
  return x
end

function RegularizedLeastSquares.normalize(::SystemMatrixBasedNormalization, op::OP, x) where {T, V <: AbstractGPUArray{T}, OP <: DenseMultiPatchOperator{T, V}}
  weights = one(real(eltype(op)))
  energy = normalize_dense_op(op, weights)
  return norm(energy)^2/size(op, 2)
end

function RegularizedLeastSquares.normalize(::SystemMatrixBasedNormalization, prod::ProdOp{T, <:WeightingOp, OP}, x) where {T, V <: AbstractGPUArray{T}, OP <: DenseMultiPatchOperator{T, V}}
  op = prod.B
  weights = prod.A.weights
  energy = normalize_dense_op(op, weights)
  return norm(energy)^2/size(prod, 2)
end

function normalize_dense_op(op::DenseMultiPatchOperator{T, V}, weights) where {T, V <: AbstractGPUArray{T}}
  backend = get_backend(op.S)
  kernel = normalize_kernel!(backend, 512)
  energy = KernelAbstractions.zeros(backend, real(eltype(op)), size(op, 1))
  kernel(energy, weights, op.S, op.xss, op.sign, Int32(div(op.M, op.nPatches)), op.RowToPatch, op.patchToSMIdx; ndrange = (512, size(op, 1)))
  return energy
end

# The normalization kernels are structured the same as the mul!-kernel. The multiplication with x is replaced by abs2 for the rownorm²
@kernel cpu = false inbounds = true function normalize_kernel!(energy, weights, @Const(S), @Const(xss), @Const(signs), @Const(M), @Const(RowToPatch), @Const(patchToSMIdx))
  # Each group/block handles a single row of the operator
  operator_row = @index(Group, Linear) # k
  patch = RowToPatch[operator_row] # p
  patch_row = mod1(operator_row, M) # j
  smIdx = patchToSMIdx[patch]
  sign = eltype(energy)(signs[patch_row, smIdx])
  @uniform grid_stride = prod(@groupsize())
  N = Int32(size(xss, 1))
  
  localIdx = @index(Local, Linear)
  shared = @localmem eltype(energy) grid_stride
  shared[localIdx] = zero(eltype(energy))

  tmp = zero(eltype(energy))
  @unroll for i = localIdx:grid_stride:N
    tmp += abs2(sign * S[xss[i, patch], patch_row, smIdx])
  end
  shared[localIdx] = tmp
  @synchronize

  @private s = div(min(grid_stride, N), Int32(2))
  while s > Int32(0)
    if localIdx <= s
      shared[localIdx] = shared[localIdx] + shared[localIdx + s]
    end
    s >>= 1
    @synchronize
  end

  if localIdx == 1
    energy[operator_row] = sqrt(get_kernel_weights(weights, operator_row)^2 * shared[localIdx])
  end
end

@inline get_kernel_weights(weights::AbstractArray, operator_row) = weights[operator_row]
@inline get_kernel_weights(weights::Number, operator_row) = weights

function Base.hash(op::DenseMultiPatchOperator{T, V}, h::UInt64) where {T, V <: AbstractGPUArray{T}}
  @warn "Hashing of GPU DenseMultiPatchOperator is inefficient"
  h = hash(typeof(op), h)
  h = @allowscalar hash(op.S, h)
  h = hash(op.grid, h)
  h = hash(op.N, h)
  h = hash(op.M, h)
  h = @allowscalar hash(op.RowToPatch, h)
  h = @allowscalar hash(op.xcc, h)
  h = @allowscalar hash(op.xss, h)
  h = @allowscalar hash(op.sign, h)
  h = hash(op.nPatches, h)
  h = @allowscalar hash(op.patchToSMIdx, h)
end