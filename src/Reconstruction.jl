export reconstruction, imToVecIm, vecImToIm
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
           images[1].data.axes...)
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
  haskey(recoParams,:emptyMeasPath) && recoParams[:emptyMeasPath]!=nothing && (recoParams[:emptyMeasPath] = MPIFile( recoParams[:emptyMeasPath] ) )

  numReco = findReco(d,study,exp,recoParams)
  if numReco > 0
    println("Reco found")
    reco = getReco(d,study,exp, numReco)
    #c = loaddata(reco.path)
    c = loadRecoDataMDF(reco.path)
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

function reconstruction(bMeas::MPIFile; kargs...)
  bSF = MPIFile(sfPath(bMeas) )
  reconstruction(bSF, bMeas; kargs...)
end

function reconstruction(filenameMeas::AbstractString; kargs...)
  bMeas = MPIFile(filenameMeas)
  reconstruction(bMeas; kargs...)
end

function reconstruction(filenameSF::AbstractString, filenameMeas::AbstractString; kargs...)
  bSF = MPIFile(filenameSF)
  bMeas = MPIFile(filenameMeas)
  reconstruction(bSF,bMeas; kargs...)
end

function reconstruction(filenameSF::AbstractString, filenameMeas::AbstractString, freq::Array; kargs...)
  bSF = MPIFile(filenameSF)
  bMeas = MPIFile(filenameMeas)
  reconstruction(bSF,bMeas,freq; kargs...)
end

function reconstruction(bSF::Union{T,Vector{T}}, bMeas::MPIFile; kargs...) where {T<:MPIFile}

  if acqNumPeriodsPerFrame(bMeas) > 1
    return reconstructionMultiPatch(bSF, bMeas; kargs...)
  else
    return reconstructionSinglePatch(bSF, bMeas;  kargs...)
  end
end

function reconstructionSinglePatch(bSF::Union{T,Vector{T}}, bMeas::MPIFile;
  minFreq=0, maxFreq=1.25e6, SNRThresh=-1,maxMixingOrder=-1, numUsedFreqs=-1, sortBySNR=false, recChannels=1:numReceivers(bMeas),
  bEmpty = nothing, bgFrames = 1, fgFrames = 1, varMeanThresh = 0, minAmplification=2, kargs...) where {T<:MPIFile}

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

function reconstruction(S, u::Array, shape; sparseTrafo = nothing,
                        lambd=0, progress=nothing, solver = "kaczmarz",
                        weights=nothing, profileName="", profiling=nothing,
                        reshapesolution = true, kargs...)

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

    if B != nothing
      d[:] = B*d #backtrafo from dual space
    end

    #if typeof(B)==LinearSolver.DSTOperator
    #	d=onGridReverse(d,shape)
    #end
    c[:,l] = real( d ) # this one is allocating
    next!(p)
    sleep(0.001)
  end

  if reshapesolution
    c = reshape(c, shape..., L)
  end
  return c
end

function reconstruction(bSF::Union{T,Vector{T}}, bMeas::MPIFile, freq::Array;
  bEmpty = nothing, bgFrames = 1,  denoiseWeight = 0, redFactor = 0.0, thresh = nothing,
  loadasreal = false, solver = "kaczmarz", sparseTrafo = nothing, saveTrafo=false,
  gridsize = gridSizeCommon(bSF), fov=calibFov(bSF), center=[0.0,0.0,0.0], useDFFoV=false,
  deadPixels=Int[], kargs...) where {T<:MPIFile}

  (typeof(bgFrames) <: AbstractRange && bEmpty==nothing) && (bEmpty = bMeas)
  bgcorrection = bEmpty != nothing ? true : false

  consistenceCheck(bSF, bMeas)

  println("Loading System matrix ...")
  S, grid = getSF(bSF, freq, sparseTrafo, solver; bgcorrection=bgcorrection, loadasreal=loadasreal,
            thresh=thresh, redFactor=redFactor, saveTrafo=saveTrafo, useDFFoV=useDFFoV,
            gridsize=gridsize, fov=fov, center=center, deadPixels=deadPixels)
  println("Loading System matrix done!")

  if denoiseWeight > 0 && sparseTrafo == nothing
    denoiseSF!(S, shape, weight=denoiseWeight)
  end

  return reconstruction(S, bSF, bMeas, freq, grid, bEmpty=bEmpty, bgFrames=bgFrames,
                        sparseTrafo=sparseTrafo, loadasreal=loadasreal,
                        solver=solver; kargs...)
end

function reconstruction(S, bSF::Union{T,Vector{T}}, bMeas::MPIFile, freq::Array, grid;
  frames = nothing, bEmpty = nothing, bgFrames = 1, nAverages = 1, sparseTrafo = nothing, loadasreal = false, maxload = 100, maskDFFOV=false,
  weightType=WeightingType.None, weightingLimit = 0, solver = "kaczmarz", spectralCleaning=true, fgFrames=1:10,
  noiseFreqThresh=0.0, kargs...) where {T<:MPIFile}

  # (typeof(bgFrames) <: AbstractRange && bEmpty==nothing) && (bEmpty = bMeas)
  bgcorrection = bEmpty != nothing ? true : false

  println("Loading emptymeas ...")
  if bEmpty!=nothing
    if acqNumBGFrames(bEmpty) > 0
      uEmpty = getMeasurementsFD(bEmpty, false, frequencies=freq, frames=measBGFrameIdx(bEmpty),
           numAverages=acqNumBGFrames(bEmpty), bgCorrection=false, loadasreal=loadasreal, spectralLeakageCorrection=spectralCleaning)
    else
      uEmpty = getMeasurementsFD(bEmpty, frequencies=freq, frames=bgFrames, numAverages=length(bgFrames),
           loadasreal=loadasreal,spectralLeakageCorrection=spectralCleaning)
    end
  end

  frames == nothing && (frames = 1:acqNumFrames(bMeas))

  weights = getWeights(weightType, freq, S, weightingLimit=weightingLimit,
                       bEmpty = bEmpty, bMeas = bMeas, bgFrames=bgFrames, bSF=bSF)

  L = -fld(-length(frames),nAverages)
  p = Progress(L, 1, "Reconstructing data...")

  #initialize output
  image = initImage(bSF,bMeas,L,grid,false)

  index = initIndex(bSF)
  iterator = nAverages == 1 ? splitrange(frames,maxload) : splitrange(frames,nAverages*maxload)
  for partframes in iterator
    println("Loading measurements ...")
    u = getMeasurementsFD(bMeas, frequencies=freq, frames=partframes, numAverages=nAverages, loadasreal=loadasreal, spectralLeakageCorrection=spectralCleaning)
    bEmpty!=nothing && (u = u .- uEmpty)

    noiseFreqThresh > 0 && setNoiseFreqToZero(u, freq, noiseFreqThresh, bEmpty = bEmpty, bMeas = bMeas, bgFrames=bgFrames)

    println("Reconstruction ...")
    c = reconstruction(S, u, shape(grid); sparseTrafo=sparseTrafo, progress=p, weights=weights,
                       reshapesolution=false, solver=solver, kargs...)
    #maskDFFOV && (c .*= trustedFOVMask(bSF))
    println("Reconstruction done")

    index = writePartToImage(c, image, index, partframes, nAverages, shape(grid))
  end

  im = vecImToIm(image)
  return im
end

function initImage(bSFFF::MPIFile, bMeas::MPIFile, L::Int, grid::RegularGridPositions,
                    loadOnlineParams=false)
  shp = shape(grid)
  T = Float32
  pixspacing = spacing(grid) ./ acqGradient(bMeas)[1] .* acqGradient(bSFFF)[1]
  offset = ffPos(bMeas) .- 0.5.*calibFov(bSFFF) .+ 0.5.*pixspacing

  Arr=Array{T}(undef, shp...,L)

  x=Axis{:x}(range(offset[1],step=pixspacing[1],length=shp[1]))
  y=Axis{:y}(range(offset[2],step=pixspacing[2],length=shp[2]))
  z=Axis{:z}(range(offset[3],step=pixspacing[3],length=shp[3]))
  t=Axis{:time}(range(0.0,step=dfCycle(bMeas),length=L))
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

function initImage(bSFs::Vector{T}, bMeas::MPIFile, L::Int, grid::RegularGridPositions,
                               loadOnlineParams=false) where {T<:MPIFile}
  return Images.ImageMeta[initImage(bSF,bMeas,L,grid,loadOnlineParams) for bSF in bSFs]
end

initIndex(bSF::MPIFile) = 1
initIndex(bSFs::Vector{T}) where {T<:MPIFile} = Int[1 for bSF in bSFs]

function writePartToImage(c, image, index::Int, partframes, nAverages, shape)
  inc = -fld(-length(partframes),nAverages)*prod(shape)
  image[index:index+inc-1] = c[:]
  index += inc
  return index
end

function writePartToImage(c, image, bSF::MPIFile)
  image[:] = c[:]
end

function writePartToImage(c, image, bSFs::Vector{T}) where {T<:MPIFile}
  #error("Where am I called")
  split = 1
  for (i,bSF) in enumerate(bSFs)
    ctemp = vec(c[split:split-1+prod(gridSize(bSF)),:])

    image[i][:] = ctemp
    split += prod(gridSize(bSF))
  end
end

function writePartToImage(c, image, index::Vector{Int}, partframes,
                     nAverages, shape)
  inc = -fld(-length(partframes),nAverages)*prod(shape)
  split = 1
  L = div(size(c,1),prod(shape))
  for i=1:L
    ctemp = c[split:split-1+prod(shape),:][:]

    image[i][index[i]:index[i]+inc-1] = ctemp
    split += prod(shape)
    index[i] += inc
  end
  return index
end

splitrange(r::Int,maxlength::Int) = r

function splitrange(r::AbstractRange,maxlength::Int)
  rout = AbstractRange[]
  stepsize = step(r)
  i = first(r)
  while i <= last(r)
    l = min(maxlength, div(last(r)-i,stepsize)+1)
    push!(rout, range(i,step=stepsize,length=l))
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
