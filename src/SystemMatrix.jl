import Base.length, Base.size

export getSF, getSparseSF, getSFGridSize, basisTrafo, deltaSampleConcentration,
       deltaSampleVolume, SVD, tikhonovLU, setlambda, CompressionAnalysis,
       idxToPos, posToIdx

using Interpolations

###############################
# delta sample functions
###############################

function deltaSampleConcentration(b::BrukerFile)
  tmp = b["PVM_MPI_TracerConcentration"]
  if tmp != nothing
    return parse(Float64, tmp)
  else
    return 1.0
  end
end

deltaSampleConcentration{T<:BrukerFile}(b::Array{T,1}) = map(deltaSampleConcentration, b)

function deltaSampleVolume(b::BrukerFile)
  V = parse(Float64, b["PVM_MPI_TracerVolume"] )*1e-6 # mu l
  return V
end

function converttoreal{T}(S::AbstractArray{Complex{T}},f)
  N = prod(calibSize(f))
  M = div(length(S),N)
  S = reshape(S,N,M)
  S = reinterpret(T,S,(2*N,M))
  p = Progress(M, 1, "Converting system matrix to real...")
  for l=1:M
    tmp = S[:,l]
    S[1:N,l] = tmp[1:2:end]
    S[N+1:end,l] = tmp[2:2:end]
    next!(p)
  end
  return S
end

function converttoreal{T}(S::AbstractArray{Complex{T},2})
  N = size(S,1)
  M = size(S,2)
  S = reinterpret(T,S,(2*N,M))
  p = Progress(M, 1, "Converting system matrix to real...")
  for l=1:M
    tmp = S[:,l]
    S[1:N,l] = tmp[1:2:end]
    S[N+1:end,l] = tmp[2:2:end]
    next!(p)
  end
  return S
end

# Systemfunction types for better inversion algorithms
@compat abstract type Systemfunction end

#type for singular value decomposition
@doc "This Type stores the singular value decomposition of a Matrix" ->
type SVD <: Systemfunction
  U::Matrix
  Σ::Vector
  V::Matrix
  D::Vector
end

SVD(U::Matrix,Σ::Vector,V::Matrix) = SVD(U,Σ,V,1./Σ)

size(A::SVD) = (size(A.U,1),size(A.V,1))
length(A::SVD) = prod(size(A))

@doc "This function can be used to calculate the singular values used for Tikhonov regularization." ->
function setlambda(A::SVD, λ::Real)
  for i=1:length(A.Σ)
    σi = A.Σ[i]
    A.D[i] = σi/(σi*σi+λ*λ)
  end
  return nothing
end

#type for Gauß elimination
type tikhonovLU <: Systemfunction
  S::Matrix
  LUfact
end

tikhonovLU(S::AbstractMatrix) = tikhonovLU(S, lufact(S'*S))

size(A::tikhonovLU) = size(A.S)
length(A::tikhonovLU) = length(A.S)

@doc "This function can be used to calculate the singular values used for Tikhonov regularization." ->
function setlambda(A::tikhonovLU, λ::Real)
  A.LUfact = lufact(A.S'*A.S + λ*speye(size(A,2),size(A,2)))
  return nothing
end

setlambda(S::AbstractMatrix, λ::Real) = nothing

function getSF(bSF, frequencies, sparseTrafo, solver; kargs...)
  if solver == "kaczmarz"
    # We now use a transpose type to simulate a row major array
    return transp(getSF(bSF, frequencies, sparseTrafo; kargs...))
  elseif solver == "pseudoinverse"
    return SVD(svd(getSF(bSF, frequencies, sparseTrafo; kargs...).')...)
  elseif solver == "cgnr" || solver == "lsqr" || solver == "fusedlasso"
    println("solver = $solver")
    return getSF(bSF, frequencies, sparseTrafo; kargs...).'
  elseif solver == "direct"
    return tikhonovLU(getSF(bSF, frequencies, sparseTrafo; kargs...).')
  else
    return getSF(bSF, frequencies, sparseTrafo; kargs...)
  end
end

# type for grid
# constructors only initializes grid
# for correct calculation call calcGrid

@compat abstract type AbstractGrid end

type KartGrid <: AbstractGrid
    Nx::Int
    Ny::Int
    Nz::Int
    X::Array{Float64,1}
    Y::Array{Float64,1}
    Z::Array{Float64,1}
end

KartGrid(g::Array{Int}) = KartGrid(g[1],g[2],g[3],0:g[1]-1,0:g[2]-1,0:g[3]-1)

function KartGrid(g::Array{Int},targetGridSize::Array{Int})
  return KartGrid(g[1],g[2],g[3],linspace(1,targetGridSize[1],g[1]),
                            linspace(1,targetGridSize[2],g[2]),
                            linspace(1,targetGridSize[3],g[3]))
end


function KartGrid(calibSize::Vector{Int},targetGridSize,
                  fov, targetFov,
                  center, targetCenter)

  minIdx = posToIdx(idxToPos([1,1,1], calibSize, fov, center),
                    targetGridSize, targetFov, targetCenter)

  maxIdx = posToIdx(idxToPos(calibSize, calibSize, fov, center),
                    targetGridSize, targetFov, targetCenter)

  for d=1:3
    if calibSize[d] == 1
      minIdx[d] = maxIdx[d] = 1
    end
  end

  grid =  KartGrid(calibSize[1],calibSize[2],calibSize[3],linspace(minIdx[1],maxIdx[1],calibSize[1]),
                            linspace(minIdx[2],maxIdx[2],calibSize[2]),
                            linspace(minIdx[3],maxIdx[3],calibSize[3]))

  return grid
end

function idxToPos(idx, size, fov, center)
  return 0.5.*fov.*(-1 + (2.*idx .- 1) ./ size) .+ center
end


function posToIdx(pos, size, fov, center)
  return 0.5*(size.* ((pos .- center) ./ ( 0.5.*fov ) + 1) + 1)
end


type ChebGrid <: AbstractGrid
    Nx::Int
    Ny::Int
    Nz::Int
    X::Array{Float64,1}
    Y::Array{Float64,1}
    Z::Array{Float64,1}
    ChebGrid(g::Array{Int})=new(g[1],g[2],g[3],0:g[1]-1,0:g[2]-1,0:g[3]-1)
end

function calcGrid(bSF::MPIFile,center::Array,voxels::Array,basisTrafo)
    pixelsDfFov=dfGrid(bSF)
    if typeof(basisTrafo)==LinearSolver.DSTOperator
        g=ChebGrid(voxels)
        g.X=cos((g.X+0.5).*pi./g.Nx)
        g.Y=cos((g.Y+0.5).*pi./g.Ny)
        g.Z=cos((g.Z+0.5).*pi./g.Nz)
    else
        g=KartGrid(voxels)
        g.X=(2.*g.X+1)./g.Nx-1
        g.Y=(2.*g.Y+1)./g.Ny-1
        g.Z=(2.*g.Z+1)./g.Nz-1
    end
    g.X=g.X.*pixelsDfFov[1]./2.+center[1]
    g.Y=g.Y.*pixelsDfFov[2]./2.+center[2]
    g.Z=g.Z.*pixelsDfFov[3]./2.+center[3]
    return g
end

# type for analysis of compression

type CompressionAnalysis
    sigma::Float64
    elements::Float64
    sigmaFreq::Array{Float64}
    CompressionAnalysis(l::Int)=new(0,0,zeros(Float64,l))
end

# this function has to be called at end of compression
function doEvaluation(compAna,sparse,nFreq,Nfull)
    if compAna!=nothing
        compAna.elements=nnz(sparse)/(nFreq*Nfull)
        compAna.sigma=mean(compAna.sigmaFreq)
    end
end

# this function has to be called in each step of compression
function calcSigmaStep(compAna,buffer,ind,n)
    if compAna!=nothing
        full=1:length(buffer)
        ind=setdiff(full,ind)
        compAna.sigmaFreq[n]=(norm(buffer[ind]))/norm(buffer)
    end
end

# calc grid of dfFov
function dfGrid(bSF::MPIFile)
    dfFovSize=dfFov(bSF)
    fovSize=calibFov(bSF)
    gridSF=calibSize(bSF)
    # adding one and floor afterwards to avoid zeros in the grid
    return floor(Int,gridSF.*dfFovSize./fovSize+1)
end
dfGrid{T<:MPIFile}(bs::Vector{T}) = [dfGrid(b) for b in bs]

# TODO TK: move this to Image Processing / Registration.jl file
# evaluates SF on the given grid
# to do the interpolation correctly a reshape on the calibSize of the systemfunction is needed (given by gridSF)
function onGrid(SF::Array,gridSF::Array{Int},grid::AbstractGrid)
  A=reshape(SF,gridSF[1],gridSF[2],gridSF[3])
  tmp = ( gridSF[1]==1 ? NoInterp() : BSpline(Linear()),
          gridSF[2]==1 ? NoInterp() : BSpline(Linear()),
          gridSF[3]==1 ? NoInterp() : BSpline(Linear()) )

  itp=interpolate(A,tmp,OnCell())

  AInterp = zeros(eltype(SF),grid.Nx,grid.Ny,grid.Nz)
  for nz=1:grid.Nz
    for ny=1:grid.Ny
      for nx=1:grid.Nx
        if grid.X[nx] >=1 && grid.X[nx] <= gridSF[1] &&
           grid.Y[ny] >=1 && grid.Y[ny] <= gridSF[2] &&
           grid.Z[nz] >=1 && grid.Z[nz] <= gridSF[3]
          AInterp[nx,ny,nz] = itp[grid.X[nx],grid.Y[ny],grid.Z[nz]]
          #if isnan(AInterp[nx,ny,nz])
          #  println("$nx $ny $nz")
          #end
        end
      end
    end
  end

  return reshape(AInterp,grid.Nx*grid.Ny*grid.Nz)
end

function onGridReverse(SF::Array,voxels)
	g=ChebGrid(voxels)
	# revert order of chebychev grid points to make them ascending
  g.X=(cos((g.X+0.5).*pi./g.Nx))[g.Nx:-1:1]
  g.Y=(cos((g.Y+0.5).*pi./g.Ny))[g.Ny:-1:1]
  g.Z=(cos((g.Z+0.5).*pi./g.Nz))[g.Nz:-1:1]
	SF=SF[g.Nx*g.Ny*g.Nz:-1:1]
	A=reshape(SF,g.Nx,g.Ny,g.Nz)
  itp=interpolate((g.X,g.Y,g.Z),A,Gridded(Linear()))
  return reshape(itp[linspace(-1,1,g.Nx),linspace(-1,1,g.Ny),linspace(-1,1,g.Nz)],g.Nx*g.Ny*g.Nz)
end

function getSFGridSize(bSF::MPIFile, sparseTrafo::Void, gridsize) #non sparse case
  return collect(gridsize)
end

function getSFGridSize(bSF::MPIFile, sparseTrafo, gridsize)
  return dfGrid(bSF)
end

function getSFGridSize{T<:MPIFile}(bSF::Vector{T}, sparseTrafo, gridsize)
  return Int64[length(bSF), getSFGridSize(bSF[1],sparseTrafo, gridsize)...]
end

function getSFVoxelSize(bSF::MPIFile, sparseTrafo::Void, gridsize) #non sparse case
  return calibFov(bSF)./collect(gridsize)
end

function getSFVoxelSize(bSF::MPIFile, sparseTrafo, gridsize)
  return voxelSize(bSF) #???
end

function getSFVoxelSize{T<:MPIFile}(bSF::Vector{T}, sparseTrafo, gridsize)
  return getSFVoxelSize(bSF[1],sparseTrafo,gridsize)
end


getSF{T<:MPIFile}(bSF::Union{T,Vector{T}}, frequencies, sparseTrafo::Void; kargs...) = getSF(bSF, frequencies; kargs...)

# TK: The following is a hack since colored and sparse trafo are currently not supported
getSF{T<:MPIFile}(bSF::Vector{T}, frequencies, sparseTrafo::AbstractString; kargs...) = getSF(bSF[1],
   frequencies, sparseTrafo; kargs...)

getSF{T<:MPIFile}(bSFs::Vector{T}, frequencies; kargs...) = cat(1,[getSF(bSF, frequencies; kargs...) for bSF in bSFs]...)

getSF(f::MPIFile; recChannels=1:numReceivers(f), kargs...) = getSF(f, filterFrequencies(f, recChannels=recChannels); kargs...)

function repairDeadPixels(S, shape, deadPixels)
  shapeT = tuple(shape...)
  #println(shapeT)
  for k=1:size(S,2)
    for dp in deadPixels
      ix,iy,iz = ind2sub(shapeT,dp)
      if 1<ix<shape[1] && 1<iy<shape[2] && 1<iz<shape[3]
         #println("$ix $iy $iz")
         lval = S[ sub2ind(shapeT,ix,iy+1,iz)  ,k]
         rval = S[ sub2ind(shapeT,ix,iy-1,iz)  ,k]
         S[dp,k] = 0.5*(lval+rval)
         fval = S[ sub2ind(shapeT,ix+1,iy,iz)  ,k]
         bval = S[ sub2ind(shapeT,ix-1,iy,iz)  ,k]
         S[dp,k] = 0.5*(fval+bval)
         uval = S[ sub2ind(shapeT,ix,iy,iz+1)  ,k]
         dval = S[ sub2ind(shapeT,ix,iy,iz-1)  ,k]
         S[dp,k] = 0.5*(uval+dval)
      end
    end
  end

end

function getSF(bSF::MPIFile, frequencies; returnasmatrix = true, procno::Integer=1,
               bgcorrection=false, loadasreal=false, loadas32bit=true,
               gridsize=calibSize(bSF), fov=calibFov(bSF), center=[0.0,0.0,0.0],
               deadPixels=Int[], kargs...)

  nFreq = rxNumFrequencies(bSF)
  shape = getSFGridSize(bSF,nothing,gridsize)

  S = getSystemMatrix(bSF, frequencies, loadas32bit=loadas32bit, bgCorrection=bgcorrection)

  if !isempty(deadPixels)
    repairDeadPixels(S,shape,deadPixels)
  end

  if collect(gridsize) != collect(calibSize(bSF)) ||
     center != [0.0,0.0,0.0] ||
     fov != calibFov(bSF)
    println("Do SF Interpolation...")
    println(gridsize, fov, center, calibSize(bSF), calibFov(bSF))

    grid = KartGrid(shape, calibSize(bSF), fov, calibFov(bSF), center, [0.0,0.0,0.0])
    #grid = KartGrid(shape, calibSize(bSF))

    SInterp = zeros(eltype(S),prod(shape),length(frequencies))
    for k=1:length(frequencies)
      SInterp[:,k] = onGrid(S[:,k],calibSize(bSF),grid)
    end
    S = SInterp
  end

  if loadasreal
    S = converttoreal(S,bSF)
    resSize = [shape..., 2*length(frequencies)]
  else
    resSize = [shape..., length(frequencies)]
  end

  if returnasmatrix
    return reshape(S,prod(resSize[1:end-1]),resSize[end])
  else
    return squeeze(reshape(S,resSize...))
  end
end



##### Sparse getSF #########

function getSF(bSF::MPIFile, frequencies, sparseTrafo::AbstractString; saveTrafo::Bool=false,
               redFactor=0.1,procno::Integer=1,  kargs...)
	if saveTrafo
		# transform whole matrix and store
		println("transformate whole matrix and store")
		basisTrafo(bSF, sparseTrafo; procno=procno, redFactor=redFactor, kargs...)
		return getSparseSF(bSF, frequencies, sparseTrafo; procno=procno, redFactor=redFactor, kargs...)
	else
		# transform only requested parts of matrix, don't store
		println("transformate only selected parts, don't store")
		return transformAndGetSparseSF(bSF, frequencies, sparseTrafo; procno=procno, redFactor=redFactor, kargs...)
	end
end

# transforms selected frequencies, do compression end return sparse matrix
# parameter globalComp: if true, each frequency gets the same number of non-zero entries
#						if false, the number of non-zero entries is determined by a changing threshold
function transformAndGetSparseSF(bSF::MPIFile,frequencies,sparseTrafo::String;voxels::Array=dfGrid(bSF),globalComp::Bool=true,redFactor=0.1,procno::Integer=1,bgcorrection=false, loadas32bit=false, loadasreal=false, compAna=nothing, combine=true, useCOM=false, depth=4, kargs ...)

    basisTrafo=linearOperator(sparseTrafo,voxels)
    z=findCenterOfDfFov(bSF;dispRes=false,combine=combine,useCOM=useCOM,depth=depth)
    grid=calcGrid(bSF,z,voxels,basisTrafo)

    localSFFilename = bgcorrection ? "systemMatrixBG" : "systemMatrix"
    sfFilename = joinpath(bSF.path,"pdata", string(procno), localSFFilename)

    nFreq = rxNumFrequencies(bSF)*rxNumChannels(bSF)
    Nfull = prod(calibSize(bSF))
    N = grid.Nx*grid.Ny*grid.Nz
    NRed = max(1,floor(Int,redFactor*N))
    l = length(frequencies)

    buffer = zeros(Complex128, N)
    indices = [Array(Int32,1) for i=1:l]
    data = [Array(Complex64,1) for i=1:l]
    numCoeff = zeros(Int32,l)
    gridSF = calibSize(bSF)

    if isfile(sfFilename)
        SF = Rawfile(sfFilename, Complex128, [Nfull,nFreq], extRaw="")
        p = Progress(nFreq, 1, "Applying basis trafo...")
        for k = 1:l
            buffer = onGrid(SF[:,frequencies[k]],gridSF,grid)
            A_mul_B!(basisTrafo,buffer)
            # compression
            if globalComp
                indices[k] = round(Int32,flipdim(sortperm(abs(buffer)),1)[1:NRed])
            else
                max, = findmax(abs(buffer))
                t = redFactor*max
                indices[k] = find(x->x>t,abs(buffer))
            end
            numCoeff[k] = length(indices[k])
            data[k] = convert(Array{Complex64,1},(buffer[indices[k]]))
            calcSigmaStep(compAna,buffer,indices[k],k)
        end
    end
    sparse = loadsparsedata(bSF,data,indices,l,N,numCoeff,loadas32bit,loadasreal)
    doEvaluation(compAna,sparse,nFreq,Nfull)
    return sparse
end

# Apply a basis transformation to the system matrix
function basisTrafo(bSF::BrukerFile, sparseTrafo::String; procno::Integer=1, redFactor=0.1, globalComp::Bool=true,voxels::Array=dfGrid(bSF),combine=true, useCOM=false, depth=4)

	path = joinpath(bSF.path,"pdata", string(procno))
	sfFilename = [joinpath(path, "systemMatrix"), joinpath(path, "systemMatrixBG")]
	sfBTFilename = [joinpath(path,"systemMatrix"*sparseTrafo),joinpath(path,"systemMatrixBG"*sparseTrafo)]
	sfIdxFilename = [joinpath(path, "systemMatrix"*sparseTrafo*"Idx"),joinpath(path, "systemMatrixBG"*sparseTrafo*"Idx")]
	sfNCFilename = [joinpath(path, "systemMatrix"*sparseTrafo*"NC"),joinpath(path, "systemMatrixBG"*sparseTrafo*"NC")]

	nFreq = rxNumFrequencies(bSF)*rxNumChannels(bSF)
	gridSF = calibSize(bSF)

	cacheAv = [ isfile(sfIdxFilename[1]) && isfile(sfBTFilename[1]) && isfile(sfNCFilename[1]),
				isfile(sfIdxFilename[2]) && isfile(sfBTFilename[2]) && isfile(sfNCFilename[2])]

	if cacheAv[1] || cacheAv[2]
		println("Basis trafo file already exists! Using the cached version.")
		return nothing
	end

	basisTrafo=linearOperator(sparseTrafo,voxels)
    z=findCenterOfDfFov(bSF;dispRes=false,combine=combine,useCOM=useCOM,depth=depth)
    grid=calcGrid(bSF,z,voxels,basisTrafo)

	Nfull = prod(gridSF)
	N=grid.Nx*grid.Ny*grid.Nz
	NRed = max(1,floor(Int,redFactor*N))

	buffer = zeros(Complex128, N)

	for l=1:2
		if isfile(sfFilename[l])
			SF = Rawfile(sfFilename[l], Complex128, [Nfull,nFreq], extRaw="")

			fdOut = open(sfBTFilename[l],"w")
			fdIndex = open(sfIdxFilename[l],"w")
			fdNC = open(sfNCFilename[l],"w")
			write(fdNC,Int32(N))

			p = Progress(nFreq, 1, "Applying basis trafo...")
			for k=1:nFreq
				buffer = onGrid(SF[:,k],gridSF,grid)
				A_mul_B!(basisTrafo, buffer)
				if globalComp
					indices = round(Int32, flipud(sortperm(abs(buffer)))[1:NRed])
				else
					max, = findmax(abs(buffer))
					t = redFactor*max
					indices = map(Int32,find(x->x>t,abs(buffer)))
				end
				d = convert(Array{Complex64,1},(buffer[indices]))
				write(fdIndex,indices)
				write(fdOut,d)
				li=length(indices)
				write(fdNC,Int32(li))
				next!(p)
			end
			close(fdOut)
			close(fdIndex)
			close(fdNC)
		end
	end
end

function loadsparsedata(f,data,indices,l,nPos,numCoeff,loadas32bit::Bool,loadasreal::Bool)
    N = nPos
    M = l
    S = loadas32bit ? map(Complex64, cat(1,data...)) : map(Complex128, cat(1,data...))
    I = round(Int64, cat(1,indices...))
    if loadasreal
        S = reshape(S,sum(numCoeff),1)
        S = converttoreal(S)
        M *= 2
        I = vcat(I,I)
        indptr = collect(0:l*2).*[0;numCoeff;numCoeff]+1
    else
        indptr = collect(0:l).*[0;numCoeff]+1
    end
    return SparseMatrixCSC(N,M,indptr,I[:],S[:])
end

getSparseSF(bSF::MPIFile, frequencies, trafoStr; redFactor=0.01, kargs...) = getSparseSF(bSF, frequencies, redFactor, trafoStr; kargs...)

# reads saved matrix from files
function getSparseSF(bSF::BrukerFile, frequencies, redFactor, trafoStr; procno::Integer=1, bgcorrection=false, loadas32bit=false, loadasreal=false, kargs ...)

	localSFFilename = bgcorrection ? "systemMatrixBG" : "systemMatrix"
	sfFilename = joinpath(bSF.path,"pdata", string(procno), localSFFilename*trafoStr)
	sfFilenameIdx = joinpath(bSF.path,"pdata", string(procno), localSFFilename*trafoStr*"Idx")
	sfFilenameNC= joinpath(bSF.path,"pdata", string(procno), localSFFilename*trafoStr*"NC")

	fs=round(Int,filesize(sfFilename)/8)
	sizeNC=round(Int,filesize(sfFilenameNC)/4)

	data = Rawfile(sfFilename, Complex64, [fs,1], extRaw="")
	indices = Rawfile(sfFilenameIdx, Int32, [fs,1], extRaw="")
	NumCoeff = Rawfile(sfFilenameNC, Int32, [sizeNC,1], extRaw="")

	NumC=NumCoeff[:,1]
	l=length(frequencies)

	i = [Array(Int32,1) for i=1:l]
	d = [Array(Complex64,1) for i=1:l]

	nPos=NumC[1,1]
	nc=reshape(NumC[2:sizeNC,1],sizeNC-1)

	for k=1:l
		if frequencies[k]==1
			lb=1
			ub=nc[frequencies[k]]
		else
			lb=sum(nc[1:(frequencies[k]-1)])+1
			ub=sum(nc[1:frequencies[k]])
		end
		i[k]=reshape(indices[lb:ub,1],ub-lb+1)
		d[k]=reshape(data[lb:ub,1],ub-lb+1)
	end

	nc=nc[frequencies]

	sparse=loadsparsedata(bSF,d,i,l,nPos,nc,loadas32bit,loadasreal)
	return sparse
end

# this function needs to be revised. Probably best if loaded from specialized file
function loadsparsedata(f,data,indices,frequencies,nPos,redFactor,nPosRed,thresh::Real,loadas32bit::Bool,loadasreal::Bool)
  warning("This function needs some work. Do not use at the Moment.")
  N = nPos
  loadasreal ? M = length(frequencies)*2 : M = length(frequencies)
  indptr = [1]
  I = Int64[]
  S = loadasreal ? loadas32bit ? Float32[] : Float64[] : loadas32bit ? Complex64[] : Complex128[]
  loadas32bit ? buffer = complex64(data[:,frequencies]) : buffer = data[:,frequencies]
  ind = indices[:,frequencies]
  p = Progress(length(frequencies), 1, "loading sparse data...")
  for k=1:length(frequencies)
    numCoeff = searchsortedlast( abs( buffer[:,k] ), map(Float32,thresh*abs( buffer[1,k] )), rev=true)
    if loadasreal
      if numCoeff != 0
        append!(S, real(buffer[1:numCoeff,k])[:])
        append!(S, imag(buffer[1:numCoeff,k])[:])
        append!(I, round(Int64,ind[1:numCoeff,k]))
        append!(I, round(Int64,ind[1:numCoeff,k]))
      end
    append!(indptr, [indptr[end]+numCoeff,indptr[end]+2*numCoeff])
    else
      if numCoeff != 0
        append!(S, buffer[1:numCoeff,k][:])
        append!(I, ind[1:numCoeff,k])
      end
      append!(indptr, [indptr[end]+numCoeff])
    end
    next!(p)
  end
  return N, M , indptr, I[:], S[:]
end
