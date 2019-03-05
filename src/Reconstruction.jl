export reconstruction
export writePartToImage, initImage

"""
This is the most high level reconstruction method using the `MDFDatasetStore`
"""
function reconstruction(d::MDFDatasetStore, study::Study, exp::Experiment, recoParams)

  !(haskey(recoParams,:SFPath)) && (recoParams[:SFPath] = sfPath( MPIFile( recoParams[:measPath] ) ))
  haskey(recoParams,:emptyMeasPath) && recoParams[:emptyMeasPath]!=nothing && (recoParams[:bEmpty] = MPIFile( recoParams[:emptyMeasPath] ) )

  numReco = findReco(d,study,exp,recoParams)
  if numReco > 0
    @info "Reconstruction found in MDF dataset store."
    reco = getReco(d,study,exp, numReco)
    c = loadRecoDataMDF(reco.path)
  else
    c = reconstruction(recoParams)
    addReco(d,study,exp, c)
  end
  return c
end

"""
This is the most high level reconstruction method that performs in-memory reconstruction
"""
function reconstruction(recoParams::Dict)
  @info "Performing in-memory reconstruction."
  bMeas = MPIFile( recoParams[:measPath] )
  !(haskey(recoParams,:SFPath)) && (recoParams[:SFPath] = sfPath( bMeas ))
  bSF = MPIFile(recoParams[:SFPath])

  c = reconstruction(bSF, bMeas; recoParams...)
  # store reco params with image
  c["recoParams"] = recoParams
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

  # Ensure that no frequencies are used that are not present in the measurement
  freq = intersect(freq, filterFrequencies(bMeas))

  @debug "selecting $(length(freq)) frequencies"

  return reconstruction(bSF, bMeas, freq; bEmpty=bEmpty, bgFrames=bgFrames, fgFrames=fgFrames, kargs...)
end


function reconstruction(bSF::Union{T,Vector{T}}, bMeas::MPIFile, freq::Array;
  bEmpty = nothing, bgFrames = 1,  denoiseWeight = 0, redFactor = 0.0, thresh = nothing,
  loadasreal = false, solver = "kaczmarz", sparseTrafo = nothing, saveTrafo=false,
  gridsize = gridSizeCommon(bSF), fov=calibFov(bSF), center=[0.0,0.0,0.0], useDFFoV=false,
  deadPixels=Int[], bgCorrectionInternal=false, kargs...) where {T<:MPIFile}

  (typeof(bgFrames) <: AbstractRange && bEmpty==nothing) && (bEmpty = bMeas)
  bgcorrection = bEmpty != nothing ? true : bgCorrectionInternal

  consistenceCheck(bSF, bMeas)

  @debug "Loading System matrix"
  S, grid = getSF(bSF, freq, sparseTrafo, solver;
            bgcorrection=bgcorrection, loadasreal=loadasreal,
            thresh=thresh, redFactor=redFactor, saveTrafo=saveTrafo,
            useDFFoV=useDFFoV, gridsize=gridsize, fov=fov, center=center,
            deadPixels=deadPixels)

  if denoiseWeight > 0 && sparseTrafo == nothing
    denoiseSF!(S, shape, weight=denoiseWeight)
  end

  return reconstruction(S, bSF, bMeas, freq, grid, bEmpty=bEmpty, bgFrames=bgFrames,
                        sparseTrafo=sparseTrafo, loadasreal=loadasreal,
                        solver=solver, bgCorrectionInternal=bgCorrectionInternal; kargs...)
end

function reconstruction(S, bSF::Union{T,Vector{T}}, bMeas::MPIFile, freq::Array, grid;
  frames = nothing, bEmpty = nothing, bgFrames = 1, nAverages = 1,
  sparseTrafo = nothing, loadasreal = false, maxload = 100, maskDFFOV=false,
  weightType=WeightingType.None, weightingLimit = 0, solver = "kaczmarz",
  spectralCleaning=true, fgFrames=1:10, bgCorrectionInternal=false,
  noiseFreqThresh=0.0, kargs...) where {T<:MPIFile}

  #(typeof(bgFrames) <: AbstractRange && bEmpty==nothing) && (bEmpty = bMeas)
  bgcorrection = bEmpty != nothing ? true : false

  @debug "Loading emptymeas ..."
  if bEmpty!=nothing
    if acqNumBGFrames(bEmpty) > 0
      uEmpty = getMeasurementsFD(bEmpty, false, frequencies=freq, frames=bgFrames, #frames=measBGFrameIdx(bEmpty),
      numAverages = acqNumBGFrames(bEmpty), bgCorrection=bgCorrectionInternal, loadasreal=loadasreal, spectralLeakageCorrection=spectralCleaning)
    else
      uEmpty = getMeasurementsFD(bEmpty, frequencies=freq, frames=bgFrames, numAverages=length(bgFrames),
      loadasreal = loadasreal,spectralLeakageCorrection=spectralCleaning)
    end
  end

  frames == nothing && (frames = 1:acqNumFrames(bMeas))

  weights = getWeights(weightType, freq, S, weightingLimit=weightingLimit,
                       bEmpty = bEmpty, bMeas = bMeas, bgFrames=bgFrames, bSF=bSF)

  L = -fld(-length(frames),nAverages) # number of tomograms to be reconstructed
  p = Progress(L, 1, "Reconstructing data...")

  #initialize output
  image = initImage(bSF,bMeas,L,nAverages,grid,false)

  currentIndex = 1
  iterator = nAverages == 1 ? Iterators.partition(frames,maxload) : Iterators.partition(frames,nAverages*maxload)
  for partframes in iterator
    @debug "Loading measurements ..."
    u = getMeasurementsFD(bMeas, frequencies=freq, frames=partframes, numAverages=nAverages, loadasreal=loadasreal, spectralLeakageCorrection=spectralCleaning, bgCorrection=bgCorrectionInternal)
    bEmpty!=nothing && (u = u .- uEmpty)

    noiseFreqThresh > 0 && setNoiseFreqToZero(u, freq, noiseFreqThresh, bEmpty = bEmpty, bMeas = bMeas, bgFrames=bgFrames)

    # convert measurement data if neccessary
    if eltype(S)!=eltype(u)
      @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
      u = map(eltype(S),u)
    end
    @debug "Reconstruction ..."
    c = reconstruction(S, u, shape(grid); sparseTrafo=sparseTrafo, progress=p,
		       weights=weights, reshapesolution=false, solver=solver, kargs...)

    currentIndex = writePartToImage!(image, c, currentIndex, partframes, nAverages)
  end

  return image
end

function writePartToImage!(image, c, currentIndex::Int, partframes, nAverages)
  # permute c's dimensions into image order
  colorsize = size(image,1)
  spatialsize = size(image,2)*size(image,3)*size(image,4)
  inc = -fld(-length(partframes),nAverages)
  c = reshape(c,spatialsize,colorsize,inc)
  c = permutedims(c,[2,1,3])
  # write c to image
  image[Axis{:time}(currentIndex:currentIndex+inc-1)] = c[:]
  currentIndex += inc
  return currentIndex
end

function initImage(bSFs::Union{T,Vector{T}}, bMeas::S, L::Int, nAverages::Int,
		   grid::RegularGridPositions, loadOnlineParams=false) where {T,S<:MPIFile}

  # the number of channels is determined by the number of system matrices
  if isa(bSFs,AbstractVector)
    numcolors = length(bSFs)
    bSF = bSFs[1]
  else
    numcolors = 1
    bSF = bSFs
  end
  # calculate axis
  shp = shape(grid)
  pixspacing = (spacing(grid) ./ acqGradient(bMeas)[1] .* acqGradient(bSF)[1])*1000u"mm"
  offset = (ffPos(bMeas) .- 0.5 .* calibFov(bSF))*1000u"mm" .+ 0.5 .* pixspacing
  dtframes = acqNumAverages(bMeas)*dfCycle(bMeas)*nAverages*1u"s"
  # initialize raw array
  Arr=Array{Float32}(undef, numcolors,shp...,L)
  # create image
  im = AxisArray(Arr, Axis{:color}(1:numcolors),
		 Axis{:x}(range(offset[1],step=pixspacing[1],length=shp[1])),
		 Axis{:y}(range(offset[2],step=pixspacing[2],length=shp[2])),
		 Axis{:z}(range(offset[3],step=pixspacing[3],length=shp[3])),
		 Axis{:time}(range(0u"ms",step=dtframes,length=L)))
  # provide meta data
  if loadOnlineParams
      imMeta = ImageMeta(im,generateHeaderDictOnline(bSF,bMeas))
  else
      imMeta = ImageMeta(im,generateHeaderDict(bSF,bMeas))
  end
  return imMeta
end

"""
Low level reconstruction method
"""
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
                            sparseTrafo=B, enforceReal=true,
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
