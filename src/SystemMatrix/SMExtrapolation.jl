export extrapolateSM

# Maybe not needed here
"""
    effrank(A)

Calculates the effective rank of a Vector or a Matrix `A`.
"""

function effrank(A)
	sig=svd(A).S
	p=sig./sum(abs.(sig))
	er=exp(-sum(p.*(log.(p))))
	return er
end

# helpful functions for building up missing data, but not essential

"""
    myellipse(N,M;a=1.0,b=1.0)

Returns all CartesianIndices of a 2D-Array of size `N`,`M` that are inside an ellipse defined by `a` and `b` around the center of the array.
"""

function myellipse(N,M;a=1.0,b=1.0)
	allInds = CartesianIndex(1,1):CartesianIndex(N,M)
	Inds = CartesianIndex[]
	mid = (1+N,1+M)./2
	for ind in allInds
		pos = ind.I .- mid
		if pos[1]^2/a^2 + pos[2]^2/b^2 <= 1
			push!(Inds,ind)
		end
	end
	return Inds
end

"""
    mycircle(N,M;rad=1.0)

Returns all CartesianIndices of a 2D-Array of size `N`,`M` that are inside a circle with radius `rad` around the center of the array.
"""

function mycircle(N,M;rad=1.0)
	return myellipse(N,M;a=rad,b=rad)
end

"""
    myellipsoid(N,M,P;a=1.0,b=1.0,c=1.0)

Returns all CartesianIndices of a 3D-Array of size `N`,`M`,`P` that are inside an ellipsoid defined by `a`, `b` and `c` around the center of the array.
"""

function myellipsoid(N,M,P;a=1.0,b=1.0,c=1.0)
	allInds = CartesianIndex(1,1,1):CartesianIndex(N,M,P)
	Inds = CartesianIndex[]
	mid = (1+N,1+M,1+P)./2
	for ind in allInds
		pos = ind.I .- mid
		if pos[1]^2/a^2 + pos[2]^2/b^2 + pos[3]^2/c^2 <= 1
			push!(Inds,ind)
		end
	end
	return Inds
end

# essential functions

"""
    get_neighbors(A,idx,neighbors)

Finds all neighbors (as defined by the input set of CartesianIndices `neighbors`) of a given set of CartesianIndices `idx` in the Array `A`.
"""

function getNeighbors(A, idx, neighbors)
	out = [i + n for i in idx for n in neighbors if i + n ∈ A]
	out = sort(unique([out; idx]))
end

"""
    getLaplacian(dims,method)

Returns the discrete laplacian in `dims=2` or `dims=3` as a (higher order block) toeplitz to perform a discrete convolution via matrix multiplication. `method=1` gives a lower order laplacian (5-point stencil in 2D, 7-point stencil in 3D), `method=2` gives a higher order laplacian.
"""

function getLaplacian(dims; method=1)
	n1,n2=dims[1:2]
	if length(dims) == 3
		n3=dims[3]
		if method == 1
			# 7-point stencil:
			D11 = spdiagm(n1,n1,-1=>ones(n1-1),0=>-2*ones(n1),1=>ones(n1-1))
			D22 = spdiagm(n2,n2,-1=>ones(n2-1),0=>-2*ones(n2),1=>ones(n2-1))
			D33 = spdiagm(n3,n3,-1=>ones(n3-1),0=>-2*ones(n3),1=>ones(n3-1))
		elseif method == 2
			# higher order stencil
			D11 = spdiagm(n1,n1,-2=>-ones(n1-2),-1=>16*ones(n1-1),0=>-30*ones(n1),
				1=>16*ones(n1-1),2=>-ones(n1-2))
			D22 = spdiagm(n2,n2,-2=>-ones(n2-2),-1=>16*ones(n2-1),0=>-30*ones(n2),
				1=>16*ones(n2-1),2=>-ones(n2-2))
			D33 = spdiagm(n3,n3,-2=>-ones(n3-2),-1=>16*ones(n3-1),0=>-30*ones(n3),
				1=>16*ones(n3-1),2=>-ones(n3-2))
		end
		return kron(I(n3),kron(I(n2),D11)) + kron(I(n3),kron(D22,I(n1))) +
				kron(kron(D33,I(n2)),I(n1))
	else
		if method == 1
			# 5-point stencil:
			D11 = spdiagm(n1,n1,-1=>ones(n1-1),0=>-2*ones(n1),1=>ones(n1-1))
			D22 = spdiagm(n2,n2,-1=>ones(n2-1),0=>-2*ones(n2),1=>ones(n2-1))
			return kron(D22,I(n1)) + kron(I(n2),D11)
		elseif method == 2
			# 9-point stencil:
			D11 = spdiagm(n1,n1,-1=>ones(n1-1),0=>10*ones(n1),1=>ones(n1-1))
			D22 = spdiagm(n2,n2,-1=>-0.5*ones(n2-1),0=>ones(n2),1=>-0.5*ones(n2-1))
			return kron(D11,D22) + kron(D22,D11)
		end
	end
end

"""
    fillmissing(A::Array;method::Integer)

Inpaints all `missing`-values in `A` using a harmonic extrapolation. `A` can be a two- or three-dimensional array. `method=1` uses a lower order laplacian (5-point stencil in 2D, 7-point stencil in 3D), `method=2` uses a higher order laplacian. An (iterative) admm-solver is used.
"""

function fillmissing(A::Array; method::Integer=1)
	# cartesian indices of unkown and known voxels
    idx_unknown = findall(@. ismissing(A))
	idx_known = findall(@. !ismissing(A))

    # build neighbor-functions 
	if ndims(A) == 2
		Ns = [[CartesianIndex((1,0)); CartesianIndex((0, 1))];
			  [CartesianIndex((1,1)); CartesianIndex((1,-1))]]
		neighbors = [Ns; -Ns]
	elseif ndims(A) == 3
		neighbors = filter(x -> x != CartesianIndex(0,0,0),
			vec([CartesianIndex(x,y,z) for x in -1:1, y in -1:1, z in -1:1]))
	end

    # get linear indices of all Neighbors of unknown voxels
	idx_work = getNeighbors(CartesianIndices(size(A)), idx_unknown, neighbors)
	lidx_work = LinearIndices(size(A))[idx_work] 

    # build Laplacian matrix
	Δ = getLaplacian(size(A); method=method)

	lidx_unknown = LinearIndices(size(A))[idx_unknown]
	lidx_known = LinearIndices(size(A))[idx_known]

    # build right hand side from known pixels
	rhs = -Δ[lidx_work, lidx_known] * A[lidx_known]

	# build solver
	reg = Regularization("L1", 0.01; shape=(length(lidx_work),length(lidx_unknown)))
	solver = createLinearSolver("admm",Δ[lidx_work, lidx_unknown];reg=reg, ρ=0.1, iterations=5)
		
	# Solving
	B = copy(A)
    # alternative IterativeSolvers.lsmr
	# B[i_unknown] .= lsmr(Δ[iwork, i_unknown], rhs; maxiter=500)
	B[lidx_unknown] .= solve(solver,convert.(Float64,rhs))

	return B
end

function extrapolateSM(bSF::MPIFile, freq::Vector{T}, ex_size; method=1, sparseTrafo=nothing, solver="kaczmarz", kargs...) where {T<:Int}
    SM, grid = getSF(bSF, freq, sparseTrafo, solver; kargs...)
    return extrapolateSM(SM, grid, ex_size; method=method)
end

function extrapolateSM(SM::AbstractMatrix, grid::RegularGridPositions, ex_size::Tuple{T,T}; method=1) where {T<:Int}
	return extrapolateSM(SM, grid, (ex_size[1],ex_size[2],0); method=method)
end

function extrapolateSM(SM::AbstractMatrix, grid::RegularGridPositions, ex_size::T; method=1) where {T<:Int}
	if shape(grid)[3] == 1
		return extrapolateSM(SM, grid, (ex_size,ex_size,0); method=method)
	else
		return extrapolateSM(SM, grid, (ex_size,ex_size,ex_size); method=method)
	end
end

function extrapolateSM(SM::AbstractMatrix, grid::RegularGridPositions, ex_size::Tuple{T,T,T}; method=1) where {T<:Int}
    
	transposed = size(SM,1) != prod(shape(grid)) ? true : false
	if transposed
		SM = transpose(SM)
	end

	K = size(SM,2)
    N1,N2,N3 = shape(grid)
    S = reshape(SM,N1,N2,N3,K)
    progress = nothing

    if N3 != 1
        M1,M2,M3 = N1+2*ex_size[1],N2+2*ex_size[2],N3+2*ex_size[3]
        S_extr = zeros(Complex{Float32},M1,M2,M3,K)
        progress==nothing ? p = ProgressMeter.Progress(K, 1, "Extrapolating 3D SystemMatrix...") : p = progress

        for k=1:K
            S_miss = convert(Array{Union{Missing, Complex{Float32}}}, S_extr[:,:,:,k])
            S_miss[2:end-1,2:end-1,2:end-1] .= missing
            S_miss[ex_size[1]+1:end-ex_size[1],ex_size[2]+1:end-ex_size[2],ex_size[3]+1:end-ex_size[3],:] = S[:,:,:,k]
            S_rl=fillmissing(real.(S_miss),method=method)
            S_im=fillmissing(imag.(S_miss),method=method)
            S_extr[:,:,:,k]=S_rl+S_im*im
            next!(p)
        end
		extrfov = (2 .* [ex_size[1], ex_size[2], ex_size[3]] .* (grid.fov ./ grid.shape)) .+ grid.fov
		extrgrid = RegularGridPositions{Float64}([M1,M2,M3], extrfov, grid.center, grid.sign)
        extrSM = transposed ? transpose(reshape(S_extr,(M1*M2*M3,K))) : reshape(S_extr,(M1*M2*M3,K))
		return extrSM,extrgrid
    else
        M1,M2 = N1+2*ex_size[1],N2+2*ex_size[2]
        S_extr = zeros(Complex{Float32},N1+2*ex_size[1],N2+2*ex_size[2],1,K)
        progress==nothing ? p = Progress(K, 1, "Extrapolating 2D SystemMatrix...") : p = progress

        for k=1:K
            S_miss = convert(Array{Union{Missing, Complex{Float32}}}, S_extr[:,:,1,k])
            S_miss[2:end-1,2:end-1] .= missing
            S_miss[ex_size[1]+1:end-ex_size[1],ex_size[2]+1:end-ex_size[2],:] = S[:,:,1,k]
            S_rl=fillmissing(real.(S_miss),method=method)
            S_im=fillmissing(imag.(S_miss),method=method)
            S_extr[:,:,1,k]=S_rl+S_im*im
            next!(p)
        end
		extrfov = (2 .* [ex_size[1], ex_size[2], 0] .* (grid.fov ./ grid.shape)) .+ grid.fov
		extrgrid = RegularGridPositions{Float64}([M1,M2,1], extrfov, grid.center, grid.sign)
		extrSM = transposed ? transpose(reshape(S_extr,(M1*M2,K))) : reshape(S_extr,(M1*M2,K))
        return extrSM,extrgrid
    end
end