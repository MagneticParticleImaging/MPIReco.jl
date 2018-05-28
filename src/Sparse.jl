export getSparseSF, getSFGridSize, basisTrafo, CompressionAnalysis


##### Sparse getSF #########

function getSF(bSF::MPIFile, frequencies, sparseTrafo::AbstractString;
               redFactor=0.1,procno::Integer=1,  kargs...)

	return transformAndGetSparseSF(bSF, frequencies, sparseTrafo; procno=procno,
     redFactor=redFactor, kargs...)
end

# transforms selected frequencies, do compression end return sparse matrix
# parameter globalComp: if true, each frequency gets the same number of non-zero entries
#						if false, the number of non-zero entries is determined by a changing threshold
function transformAndGetSparseSF(bSF::MPIFile,frequencies,sparseTrafo::String;voxels::Array=dfGrid(bSF),
  globalComp::Bool=true,redFactor=0.1,procno::Integer=1,bgcorrection=false, loadas32bit=true,
     loadasreal=false, compAna=nothing, combine=true, useCOM=false, depth=4, kargs...)

    basisTrafo = linearOperator(sparseTrafo,voxels)
    z = findCenterOfDfFov(bSF;dispRes=false,combine=combine,useCOM=useCOM,depth=depth)
    println("Found center at $z")

    grid = calcGrid(bSF,z,voxels,basisTrafo)

    nFreq = rxNumFrequencies(bSF)*rxNumChannels(bSF)
    Nfull = prod(calibSize(bSF))
    N = grid.Nx*grid.Ny*grid.Nz
    NRed = max(1,floor(Int,redFactor*N))
    l = length(frequencies)

    buffer = zeros(Complex64, N)
    indices = [zeros(Int32,1) for i=1:l]
    data = [zeros(Complex64,1) for i=1:l]
    numCoeff = zeros(Int32,l)
    gridSF = calibSize(bSF)

    p = Progress(nFreq, 1, "Applying basis trafo...")
    for k = 1:l
        SF = map(Complex64, systemMatrix(bSF, frequencies[k], bgcorrection) )
        buffer = onGrid(SF,gridSF,grid)
        A_mul_B!(basisTrafo,buffer)
        # compression
        if globalComp
            indices[k] = round.(Int32,flipdim(sortperm(abs.(buffer)),1)[1:NRed])
        else
            max, = findmax(abs.(buffer))
            t = redFactor*max
            indices[k] = find(x->x>t,abs.(buffer))
        end
        numCoeff[k] = length(indices[k])
        data[k] = map(Complex64,(buffer[indices[k]]))
        calcSigmaStep(compAna,buffer,indices[k],k)
    end

    sparse = loadsparsedata(bSF,data,indices,l,N,numCoeff,loadas32bit,loadasreal)
    doEvaluation(compAna,sparse,nFreq,Nfull)
    return sparse
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
