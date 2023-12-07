export smRecovery

###################################################################
# System Matrix recovery using DCT-sparsity and low-rank constraing
###################################################################
function smRecovery(y::Matrix{T}, samplingIdx::Array{Int64}, params::Dict) where T
  shape = params[:shape]

  # sampling operator
  P = SamplingOp(Complex{real(T)}, pattern=samplingIdx, shape=shape)

  # sparsifying transformation and regularization
  if !haskey(params, :sparseTrafo)
    sparseTrafo = DCTOp(Complex{real(T)},shape=params[:shape],dcttype=2)
  else
    sparseTrafo = params[:sparseTrafo]
  end
  reg1 = TransformedRegularization( L1Regularization(params[:λ_l1]), sparseTrafo)

  # setup low rank regularization
  λ_lr = get(params, :λ_lr, 0.0)
  ρ_lr = get(params, :ρ_lr, 0.0)
  if λ_lr==0 || ρ_lr==0
    # no LR regularization
    reg = [reg1]
    ρ = [params[:ρ_l1]]
  elseif length(shape)==2
    # 2d case (use nuclear norm)
    regNN = regularizationNN(λ_lr, shape)
    reg = [reg1, regNN]
    ρ = [params[:ρ_l1], ρ_lr]
  else
    # 3d case (use tensor nuclear norm)
    regTNN1 = regularizationTNN(λ_lr, shape, 1)
    regTNN2 = regularizationTNN(λ_lr, shape, 2)
    regTNN3 = regularizationTNN(λ_lr, shape, 3)
    reg = [reg1, regTNN1, regTNN2, regTNN3]
    ρ = [params[:ρ_l1], ρ_lr, ρ_lr, ρ_lr]
  end

  # precalculate inverse for split Bregman iteration
  precon = invPrecon(samplingIdx, prod(shape); params...)
  params[:precon] = precon
  params[:iterationsCG] = 1

  # normalized measurement
  y2, y_norm = getNormalizedMeasurement(y)

  solver = [SplitBregman(P; reg=deepcopy(reg), ρ=ρ, params...) for i=1:Threads.nthreads()]

  sfMat = zeros(ComplexF64,prod(shape),size(y,2))
  @time Threads.@threads for k=1:size(y,2)
    t = Threads.threadid()
    sfMat[:,k] .= solve(solver[t], y2[k])
    # undo normalization
    sfMat[:,k] *= y_norm[k]
  end

  return sfMat
end

#########################
# low rank regularization
#########################
function regularizationTNN(λ::AbstractFloat, shape::NTuple{3,Int64}, mode::Int64)
  proxMap = (x,y;kargs...)->proxTNN!(x,y,shape,mode)
  reg = Regularization(prox! = proxMap, λ=λ, params=Dict{Symbol,Any}())
end

function proxTNN!(x, λ::Real, shape::NTuple{3,Int64},mode::Int64; kargs...)
  nx,ny,nz = shape
  # flattened tensor
  x = reshape(x,shape)
  if mode==1
    x_flat = reshape(x,nx,ny*nz)
  elseif mode==2
    x_flat = reshape(permutedims(x,[2,3,1]),ny,nx*nz)
  else
    if mode!=3
      @error "mode in proxTNN! should be 1,2 or 3"
    end
    x_flat = reshape(permutedims(x,[3,1,2]),nz,nx*ny)
  end
  # perform singular value thresholding
  usv = svd!(x_flat)
  proxL1!(usv.S,λ)
  x_flat[:] .= vec(usv.U*Diagonal(usv.S)*usv.Vt)
  # undo flattening operation
  if mode==1
    x[:] .= vec(x_flat)
  elseif mode==2
    x[:] .= vec(permutedims(reshape(x_flat,ny,nz,nx),[3,1,2]))
  else
    x[:] .= vec(permutedims(reshape(x_flat,nz,nx,ny),[2,3,1]))
  end
end

# low rank regularization for 2d (using SVD)
function regularizationNN(λ::Real, shape::NTuple{2,Int64} )
  proxMap = (x,y;kargs...)->proxNN!(x,y,shape;kargs...)
  return Regularization(prox! = proxMap, λ=λ)
end

# singular-value thresholding
function proxNN!(x, λ::Real, shape::NTuple{2,Int64}; kargs...)
  U,S,V = svd(reshape(x,shape))
  proxL1!(S,λ)
  x[:] .= vec(U*Matrix(Diagonal(S))*V')
end

# fixed rank constraing for 3d (using HOSVD)
function FRTruncation!(x::Array{T,3}, rank::NTuple{3,Int64}) where T
  guvw_r = hosvd(real.(x),rank)
  guvw_i = hosvd(imag.(x),rank)
  x[:] .= compose(guvw_r)+1im*compose(guvw_i)
end

# fixed rank constraing for 2d (using SVD)
function FRTruncation!(x, rank::Int64)
  U,S,V = svd(reshape(x, shape))
  S[rank+1:end] .= 0.0
  x[:] .= vec(U*Matrix(Diagonal(S))*V')
end

##########################
# inverse of normal matrix
##########################
mutable struct diagPrecon{T}
  diag::Vector{T}
end

function ldiv!(y::Vector{T},P::diagPrecon{T},x::Vector{T}) where T
  y[:] .= P.diag .* x
end

function ldiv!(P::diagPrecon{T},x::Vector{T}) where T
  x[:] .= P.diag .* x
end

\(P::diagPrecon{T}, x::Vector{T}) where T = P.diag .* x

function invPrecon(samplingIdx, ncol; ρ::Vector{Float64}=[1.0], kargs...)
  normalP = normalPDiag(samplingIdx, ncol)
  diag = [normalP[i]+sum(ρ) for i=1:length(normalP)]
  diag = diag.^-1
  return diagPrecon(diag)
end

function normalPDiag(samplingIdx, ncol; T=ComplexF64)
  normalP = zeros(T,ncol)
  for i=1:length(samplingIdx)
    normalP[samplingIdx[i]] = 1.0
  end
  return normalP
end

function getNormalizedMeasurement(y::Matrix{T}) where T
  y_norm = [norm(y[:,k]) for k=1:size(y,2)]
  y2 = [ComplexF64.(y[:,k]/y_norm[k]) for k=1:size(y,2)]
  return y2, y_norm
end
