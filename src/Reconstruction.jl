export reconstruction, reconstructionSFCutOff, imToVecIm, vecImToIm, getshape
export writePartToImage, initImage

function imToVecIm(image::ImageMeta)
  out = ImageMeta[]
  for i=1:size(image,1)
    I = getindex(image, i, ntuple(x->:,ndims(image)-1)...)
    push!(out, I)
  end
  return out
end

function vecImToIm(images::Vector)
  out = AxisArray(zeros(length(images),size(images[1])...),
           AxisArrays.Axis{:color}(1:length(images)),
           Images.axes(images[1].data)...)
  for i in eachindex(images)
    out[i,ntuple(x->:,ndims(images[1]))...] = images[i]
  end
  out = ImageMeta(out,properties(images[1]))
  return out
end

function vecImToIm(image::ImageMeta)
  out = AxisArray(zeros(1, size(image)...),
           AxisArrays.Axis{:color}(1:1),
           Images.axes(image.data)...)

  out[1,ntuple(x->:,ndims(image))...] = image

  out = ImageMeta(out,properties(image))
  return out
end


function reconstruction(d::MDFDatasetStore, study::Study, exp::Experiment, recoParams)

  !(haskey(recoParams,:SFPath)) && (recoParams[:SFPath] = sfPath( MPIFile( recoParams[:measPath] ) ))

  numReco = findReco(d,study,exp,recoParams)
  if numReco > 0
    println("Reco found")
    reco = getReco(d,study,exp, numReco)
    c = loaddata(reco.path)
  else
    println("Reco not found")
    c = reconstruction(recoParams)
    addReco(d,study,exp, c)
  end
  return c
end

"""
This is the most high level reconstruction method that performs in-memory reconstruction
"""
function reconstruction(recoParams::Dict)
  bMeas = MPIFile( recoParams[:measPath] )
  !(haskey(recoParams,:SFPath)) && (recoParams[:SFPath] = sfPath( bMeas ))
  bSF = MPIFile(recoParams[:SFPath])

  c = reconstruction(bSF, bMeas; recoParams...)
  #make me collored??
  #c[1]["recoParams"] = recoParams
  if ndims(c)==1
    c[1]["recoParams"] = recoParams
  else
    c["recoParams"] = recoParams
  end

  return c
end




@doc """
This file contains MPI reconstruction algorithms. Different forms of the
reconstruction function are available.

**Optional keyword arguments**

* `lambd::Float64`: The regularization parameter, relative to the matrix trace.
* `iterations::Int`: Number of iterations of the iterative solver.
* `solver::AbstractString`: Algorithm used to solve the imaging equation (currently "kaczmarz" or "cgnr").
* `normWeights::Bool`: Enable row normalization (true/false).
* `denoiseWeight::Float`: Denoise the system matrix prior to reconstrution. Disabled if denoiseWeight=0.
* `sparseTrafo::AbstractString/Nothing`: Enable sparseTrafo if set to "DCT" or "FFT". Disable sparseTrafo if set to nothing.
""" ->
function reconstruction(bMeas::MPIFile; kargs...)
  bSF = MPIFile(sfPath(bMeas) )
  reconstruction(bSF, bMeas; kargs...)
end

@doc """
Returns a MPI image.

**Arguments**

* `filenameMeas::AbstractString´: file path to measurement.

**Optional keyword arguments**

* `kargs`: Additional arguments to be passed to low level reconstruction.

**Returns**

* `Array`: Returns Array or stack containing Arrays depending on the number of reconstructed frames.
"""->
function reconstruction(filenameMeas::AbstractString; kargs...)
  bMeas = MPIFile(filenameMeas)
  reconstruction(bMeas; kargs...)
end

@doc """
Returns a MPI image.

**Arguments**

* `filenameSF::AbstractString´: file path to system matrix.
* `filenameMeas::AbstractString´: file path to measurement.

**Optional keyword arguments**

* `kargs`: Additional arguments to be passed to low level reconstruction.

**Returns**

* `Array`: Returns Array or stack containing Arrays depending on the number of reconstructed frames.
"""->
function reconstruction(filenameSF::AbstractString, filenameMeas::AbstractString; kargs...)
  bSF = MPIFile(filenameSF)
  bMeas = MPIFile(filenameMeas)
  reconstruction(bSF,bMeas; kargs...)
end

@doc """
Returns a MPI image.

**Arguments**

* `filenameSF::AbstractString´: file path to system matrix.
* `filenameMeas::AbstractString´: file path to measurement.
* `frequencies::Array´: indices of frequencies to be used for reconstruction.

**Optional keyword arguments**

* `kargs`: Additional arguments to be passed to low level reconstruction.

**Returns**

* `Array`: Returns Array or stack containing Arrays depending on the number of reconstructed frames.
"""->
function reconstruction(filenameSF::AbstractString, filenameMeas::AbstractString, freq::Array; kargs...)
  bSF = MPIFile(filenameSF)
  bMeas = MPIFile(filenameMeas)
  reconstruction(bSF,bMeas,freq; kargs...)
end

@doc """
Returns a MPI image.

**Arguments**

* `S::AbstractArray´: Array containing system matrix.
* `u::Array´: Array containing measurements.
* `shape`: Array/Tuple containing the grid size of the system matrix.

**Optional keyword arguments**

* `weightingLimit::Float64`: Limit the range of the weights.
* `solver::AbstractString`: Algorithm used to solve the imaging equation (currently "kaczmarz" or "cgnr").
* `lambd::Float64`: The regularization parameter, relative to the matrix trace.
* `progress::Progress`: If progress is provided it is increased by one after each reconstructed image.
* `kargs`: Additional arguments to be passed to low level reconstruction.

**Returns**

* `Array`: Returns Array or stack containing Arrays depending on the number of reconstructed frames.
"""->
function reconstruction(S, u::Array, shape; sparseTrafo = nothing,
                                            lambd=0, progress=nothing, solver = "kaczmarz", weights=nothing, profileName="", profiling=nothing, reshapesolution = true, kargs...)

  N = size(S,2) #prod(shape)
  M = div(length(S), N)

  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N,L)
  #c = zeros(real(eltype(u)),N,L) Change by J.Dora

  if lambd > 0
    trace = calculateTraceOfNormalMatrix(S,weights)
    lambd *= trace / N
    setlambda(S,lambd)
  end

  #solv = linearSolver(solver)
  B = linearOperator(sparseTrafo, shape)
  solv = createLinearSolver(solver, S; shape=shape, weights=weights, lambdL2=lambd,
                            sparseTrafo=B, verbose = false, enforceReal=true,
                            enforcePositive=true, kargs...)

  progress==nothing ? p = Progress(L, 1, "Reconstructing data...") : p = progress
  for l=1:L

    d = solve(solv, u[:,l])

    A_mul_B!(B, d) #backtrafo from dual space
    if typeof(B)==LinearSolver.DSTOperator
    	d=onGridReverse(d,shape)
    end
    c[:,l] = real( d ) # this one is allocating
    next!(p)
    sleep(0.001)
  end

  if reshapesolution
    c = reshape(c, shape..., L)
  end
  return c
end


@doc """
Returns a MPI image.

**Arguments**

* `bSF::BrukerFile`: Handle to first system matrix.
* `bMeas::BrukerFile`: Handle to measurement.
* `frequencies::Array`: indices of frequencies to be used for reconstruction.

**Optional keyword arguments**

* `frames::Range{Int64}`: Index set of frames to be reconstructed .
* `bEmpty::BrukerFile`: Handle of background measurement to be subtracted from measurement prior to reconstruction.
* `denoiseWeight::Float64`: ???
* `bgFrames::Union(Range{Int64},Int64)`: Index set of frames used for background correction.
* `nAverages::Int`: ???
* `sparseTrafo::AbstractString`: Sparse transformation (e.g. DCT, FFT) used for matrix compression.
* `redFactor::Real`: Compression rate of sparse transformation.
* `thresh::Real`: Compression cuts off elements based on there realtive energy. Disables `redFactor`!
* `loadasreal:: Bool`: Load system matrix as real matrix.
* `loadas32bit:: Bool` Load system matrix with single precicion float.
* `maxload::Int`: delimits the number of measurements loaded into memory at a time.
* `kargs`: Additional arguments to be passed to low level reconstruction.

**Returns**

* `Array`: Returns Array or stack containing Arrays depending on the number of reconstructed frames.
"""->
function reconstruction{T<:MPIFile}(bSF::Union{T,Vector{T}}, bMeas::MPIFile, freq::Array;
  bEmpty = nothing, bgFrames = 1,  denoiseWeight = 0, redFactor = 0.0, thresh = nothing,
  loadasreal = false, loadas32bit = true, solver = "kaczmarz", sparseTrafo = nothing, saveTrafo=false,
  gridsize = gridSizeCommon(bSF), fov=calibFov(bSF), center=[0.0,0.0,0.0],
  deadPixels=Int[], kargs...)

  shape = getSFGridSize(bSF, sparseTrafo, gridsize)
  (typeof(bgFrames) <: Range && bEmpty==nothing) && (bEmpty = bMeas)
  bgcorrection = bEmpty != nothing ? true : false

  consistenceCheck(bSF, bMeas)

  S = getSF(bSF, freq, sparseTrafo, solver; bgcorrection=bgcorrection, loadasreal=loadasreal,
            thresh=thresh, redFactor=redFactor, saveTrafo=saveTrafo,
            gridsize=gridsize, fov=fov, center=center, deadPixels=deadPixels)

  if denoiseWeight > 0 && sparseTrafo == nothing
    denoiseSF!(S, shape, weight=denoiseWeight)
  end

  return reconstruction(S, bSF, bMeas, freq, bEmpty=bEmpty, bgFrames=bgFrames,
                        sparseTrafo=sparseTrafo, loadasreal=loadasreal,
                        solver=solver, gridsize=gridsize; kargs...)
end

function reconstruction{T<:MPIFile}(S, bSF::Union{T,Vector{T}}, bMeas::MPIFile, freq::Array;
  frames = nothing, bEmpty = nothing, bgFrames = 1, nAverages = 1, sparseTrafo = nothing, loadasreal = false, maxload = 100, maskDFFOV=false,
  weightType=WeightingType.None, weightingLimit = 0, solver = "kaczmarz", spectralCleaning=true, fgFrames=1:10,
  noiseFreqThresh=0.0, gridsize = gridSizeCommon(bSF), kargs...)

  shape = getSFGridSize(bSF, sparseTrafo, gridsize)
  # (typeof(bgFrames) <: Range && bEmpty==nothing) && (bEmpty = bMeas)
  bgcorrection = bEmpty != nothing ? true : false

  bEmpty!=nothing && (uEmpty = getMeasurementsFD(bEmpty, frequencies=freq, frames=bgFrames, numAverages=length(bgFrames),
                                                 loadasreal=loadasreal,spectralLeakageCorrection=spectralCleaning))

  frames == nothing && (frames = 1:acqNumFrames(bMeas))

  weights = getWeights(weightType, freq, S, weightingLimit=weightingLimit,
                       bEmpty = bEmpty, bgFrames=bgFrames, bMeas = bMeas, bgFrames=bgFrames, bSF=bSF)

  L = -fld(-length(frames),nAverages)
  p = Progress(L, 1, "Reconstructing data...")

  #initialize output
  image = initImage(bSF,bMeas,L,true,sparseTrafo,false,gridsize)

  index = initIndex(bSF)
  iterator = nAverages == 1 ? splitrange(frames,maxload) : splitrange(frames,nAverages*maxload)
  for partframes in iterator
    u = getMeasurementsFD(bMeas, frequencies=freq, frames=partframes, numAverages=nAverages, loadasreal=loadasreal, spectralLeakageCorrection=spectralCleaning)
    bEmpty!=nothing && (u = u .- uEmpty)

    noiseFreqThresh > 0 && setNoiseFreqToZero(u, freq, noiseFreqThresh, bEmpty = bEmpty, bgFrames=bgFrames, bMeas = bMeas, bgFrames=bgFrames)

    c = reconstruction(S, u, shape; sparseTrafo=sparseTrafo, progress=p, weights=weights,
                       reshapesolution=false, solver=solver, kargs...)
    #maskDFFOV && (c .*= trustedFOVMask(bSF))
    #c = reshape(c,prod(shape),size(c)[end])
    if size(c,1) != prod(shape)
      println("This code does not look correct")
      c = c[1:prod(shape),:]
    end
    index = writePartToImage(c, image, index, partframes, nAverages, bSF, sparseTrafo, gridsize)
  end

  im = vecImToIm(image)
  return im
end

function getImageGridSize(bSF, sparseTrafo)
  return getshape( sparseTrafo == nothing ? gridSize(bSF) : dfGrid(bSF))
end

getshape(shapes::Array{Int64,1}) = shapes#[ find(shapes .> 1)]
getshape(shapes::Array) = [length(shapes), shapes[1]...] #This was wrong: sum([prod(s) for s in shapes])

function initImage(bSFFF::MPIFile, bMeas::MPIFile, L::Int, loadas32bit::Bool,
                   sparseTrafo, loadOnlineParams=false, gridsize=gridSizeCommon(bSFFF))
  shape = getSFGridSize(bSFFF, sparseTrafo, gridsize)
  T = Float32
  pixspacing = getSFVoxelSize(bSFFF, sparseTrafo, gridsize) ./ acqGradient(bMeas)[1] .* acqGradient(bSFFF)[1]
  offset = ffPos(bMeas) .- 0.5.*calibFov(bSFFF) .+ 0.5.*pixspacing

  Arr=Array{T}(shape...,L)

  x=Axis{:x}(range(offset[1],pixspacing[1],shape[1]))
  y=Axis{:y}(range(offset[2],pixspacing[2],shape[2]))
  z=Axis{:z}(range(offset[3],pixspacing[3],shape[3]))
  t=Axis{:time}(range(0.0,dfCycle(bMeas),L))
  im = AxisArray(Arr,x,y,z,t)

  # The following does for some reason not work: ReadOnlyMemoryError
  #im = AxisArray(Arr, (:x,:y,:z,:time),tuple(pixspacing...,dfCycle(bMeas)),tuple(offset...,0.0))
  if loadOnlineParams
      imMeta = ImageMeta(im,generateHeaderDictOnline(bSFFF,bMeas))
  else
      imMeta = ImageMeta(im,generateHeaderDict(bSFFF,bMeas))
  end
  return imMeta
end

function initImage{T<:MPIFile}(bSFs::Vector{T}, bMeas::MPIFile, L::Int,
                               loadas32bit::Bool,sparseTrafo,loadOnlineParams=false, gridsize=gridSize(bSFs[1]))
  return Images.ImageMeta[initImage(bSF,bMeas,L,loadas32bit,sparseTrafo,loadOnlineParams, gridsize) for bSF in bSFs]
end

initIndex(bSF::MPIFile) = 1
initIndex{T<:MPIFile}(bSFs::Vector{T}) = Int[1 for bSF in bSFs]

function writePartToImage(c, image, index::Int, partframes, nAverages,
                          bSF::MPIFile, sparseTrafo, gridsize)
  #if deltaSampleConcentration(bSF) > 0  # this is usually not helpful
  #  scale!(c,1/deltaSampleConcentration(bSF))
  #end
  shape = getSFGridSize(bSF, sparseTrafo, gridsize)
  inc = -fld(-length(partframes),nAverages)*prod(shape)
  image[index:index+inc-1] = c[:]
  index += inc
  return index
end

function writePartToImage(c, image, bSF::MPIFile)
  #if deltaSampleConcentration(bSF) > 0
  #  scale!(c,1/deltaSampleConcentration(bSF))
  #end
  image[:] = c[:]
end


function writePartToImage{T<:MPIFile}(c, image, bSFs::Vector{T})
  #error("Where am I called")
  split = 1
  for (i,bSF) in enumerate(bSFs)
    ctemp = vec(c[split:split-1+prod(gridSize(bSF)),:])
    #if deltaSampleConcentration(bSF) > 0
    #  scale!(ctemp,1/deltaSampleConcentration(bSF))
    #end
    image[i][:] = ctemp
    split += prod(gridSize(bSF))
  end
end

function writePartToImage{T<:MPIFile}(c, image, index::Vector{Int}, partframes,
                     nAverages, bSFs::Vector{T}, sparseTrafo, gridsize)
  shape = getSFGridSize(bSFs[1], sparseTrafo, gridsize)
  inc = -fld(-length(partframes),nAverages)*prod(shape)
  split = 1
  for (i,bSF) in enumerate(bSFs)
    ctemp = c[split:split-1+prod(shape),:][:]
    #if deltaSampleConcentration(bSF) > 0
    #  scale!(ctemp,1/deltaSampleConcentration(bSF))
    #end
    image[i][index[i]:index[i]+inc-1] = ctemp
    split += prod(shape)
    index[i] += inc
  end
  return index
end

splitrange(r::Int,maxlength::Int) = r

function splitrange(r::Range,maxlength::Int)
  rout = Range[]
  stepsize = step(r)
  i = start(r)
  while i<=last(r)
    l = min(maxlength, div(last(r)-i,stepsize)+1)
    push!(rout, range(i,stepsize,l))
    i += stepsize*l
  end
  return rout
end

function splitrange(r::Array{Int,1}, maxlength::Int)
  res=Array[]
  l=length(r)
  if maxlength < l
    i=1
    while i <= l
      if (i+maxlength-1) > l
        push!(res, r[i:end])
      else
        push!(res, r[i:i+maxlength-1])
      end
      i=i+maxlength
    end
    return res
  else
    return Any[r]
  end
end

@doc """
Returns a MPI image.

**Arguments**

* `bSF::BrukerFile´: Handle to system matrix.
* `bMeas::BrukerFile´: Handle to measurement.

**Optional keyword arguments**

* `minFreq::Real`: Lower frequency used for cut off.
* `maxFreq::Real`: Higher frequency used for cut off.
* `SNRThresh::Real`: Signal to noise threshold used to filter out frequency components.
* `sortBySNR::Bool`: ???
* `recChannels::Range`: Specify recieving channels to use for reconstruction.
* `kargs`: Additional arguments to be passed to lower level reconstruction.

**Returns**

* `Array`: Returns Array or stack containing Arrays depending on the number of reconstructed frames.
"""->
function reconstructionSinglePatch{T<:MPIFile}(bSF::Union{T,Vector{T}}, bMeas::MPIFile;
  minFreq=0, maxFreq=1.25e6, SNRThresh=-1,maxMixingOrder=-1, numUsedFreqs=-1, sortBySNR=false, recChannels=1:numReceivers(bMeas),
  bEmpty = nothing, bgFrames = 1, fgFrames = 1, varMeanThresh = 0, minAmplification=2, kargs...)

  freq = filterFrequencies(bSF,minFreq=minFreq, maxFreq=maxFreq,recChannels=recChannels, SNRThresh=SNRThresh, numUsedFreqs=numUsedFreqs, sortBySNR=sortBySNR)

  if varMeanThresh > 0
    bEmptyTmp = (bEmpty == nothing) ? bMeas : bEmpty

    freqVarMean = filterFrequenciesVarMean(bMeas, bEmptyTmp, fgFrames, bgFrames;
                                        thresh=varMeanThresh, minAmplification=minAmplification,
                                        minFreq=minFreq, maxFreq=maxFreq,recChannels=recChannels)
    freq = intersect(freq, freqVarMean)
  end

  println("Frequency Selection: ", length(freq), " frequencies")

  # Ensure that no frequencies are used that are not present in the measurement
  freq = intersect(freq, filterFrequencies(bMeas))

  println("Frequency Selection: ", length(freq), " frequencies")

  return reconstruction(bSF, bMeas, freq; bEmpty=bEmpty, bgFrames=bgFrames, fgFrames=fgFrames, kargs...)
end

function reconstruction{T<:MPIFile}(bSF::Union{T,Vector{T}}, bMeas::MPIFile; kargs...)

  if acqNumPeriodsPerFrame(bMeas) > 1
    return reconstructionMultiPatch(bSF, bMeas; kargs...)
  else
    return reconstructionSinglePatch(bSF, bMeas;  kargs...)
  end
end
