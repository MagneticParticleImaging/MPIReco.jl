export reconstruction
export writePartToImage!, initImage

function reconstruction(experiment::Experiment{MDFDatasetStore}, recoParams)
    study = getStudy(experiment)
    dataSetStore = getMDFStore(experiment)
    return reconstruction(dataSetStore, study, experiment, recoParams)
end

function reconstruction(dataSetStore::MDFDatasetStore,experiment::Experiment, recoParams)
    study = getStudy(experiment)
    return reconstruction(dataSetStore, study, experiment, recoParams)
end

"""
This is the most high level reconstruction method using the `MDFDatasetStore`
"""
function reconstruction(d::MDFDatasetStore, study::Study, experiment::Experiment, recoParams)
  recoParams = copy(recoParams)
  if haskey(recoParams,:measPath)
      @error "conflicting measurement paths in `recoParams` and `experiment`. Set `recoParams[:measPath] = path(experiment)`"
  end

  recoParams[:measPath] = path(experiment)

  !(haskey(recoParams,:SFPath)) && (recoParams[:SFPath] = sfPath( MPIFile( recoParams[:measPath] ) ))
  if isa(recoParams[:SFPath], AbstractString)
    recoParams[:SFPath] =  MPIFiles.extendPath(d,recoParams[:SFPath])
  else
    recoParams[:SFPath] = map(s->MPIFiles.extendPath(d,s),recoParams[:SFPath]) 
  end
  
  numReco = findReco(d,study,experiment,recoParams)

  haskey(recoParams,:emptyMeasPath) && recoParams[:emptyMeasPath]!=nothing && (recoParams[:emptyMeas] = MPIFile( recoParams[:emptyMeasPath] ) )

  if numReco > 0
    @info "Reconstruction found in MDF dataset store."
    reco = getReco(d,study,experiment, numReco)
    c = loadRecoData(reco.path)
  else
    c = reconstruction(recoParams)
    addReco(d,study,experiment, c)
  end
  return c
end

# The previous function is somewhat redundant in its arguments. In particular
# the measPath and the study/experiment are redundant. Should be cleaned up before
# it can be used as a user facing API
#
#function reconstruction(d::MDFDatasetStore, recoParams::Dict)
#
#  study = ???
#
#  studies = getStudies( activeDatasetStore(m) )
#
#  for study in studies
#    push!(m.studyStore, (study.date, study.name, study.subject, study.path, true))
#  end
#
#  m.currentStudy = Study(TreeModel(m.studyStoreSorted)[currentIt,4],
#                         TreeModel(m.studyStoreSorted)[currentIt,2],
#                         TreeModel(m.studyStoreSorted)[currentIt,3],
#                         TreeModel(m.studyStoreSorted)[currentIt,1])
#
#  experiment = getExperiment(study, recoParams[:measPath])
#  return reconstruction(d, study, experiment, recoParams)
#end


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
  c.recoParams = recoParams
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
  if haskey(kargs, :periodicMotionCorrection) && kargs[:periodicMotionCorrection]
    return reconstructionPeriodicMotion(bSF, bMeas; kargs...)
  elseif acqNumPeriodsPerFrame(bMeas) > 1 &&
      rxNumSamplingPoints(bSF) == rxNumSamplingPoints(bMeas) &&
      (acqNumPeriodsPerFrame(bSF) == 1 || typeof(bSF) == MultiMPIFile)
    # This branch is only used of the measurements are multi-patch and if the
    # system matrix is not fully sampled. This is the case if either we have a
    # single system matrix that is reused for-multiple patches, or if we have
    # a MultiMPIFile
    return reconstructionMultiPatch(bSF, bMeas; kargs...)
  else
    return reconstructionSinglePatch(bSF, bMeas;  kargs...)
  end
end

function reconstructionSinglePatch(bSF::Union{T,Vector{T}}, bMeas::MPIFile;
  minFreq=0, maxFreq=1.25e6, SNRThresh=-1,maxMixingOrder=-1, numUsedFreqs=-1, sortBySNR=false, recChannels=1:numReceivers(bMeas),
  bEmpty = nothing, emptyMeas=bEmpty, bgFrames = 1, fgFrames = 1, varMeanThresh = 0, minAmplification=2, 
  numPeriodAverages=1, numPeriodGrouping=1, kargs...) where {T<:MPIFile}

  freq = filterFrequencies(bSF,minFreq=minFreq, maxFreq=maxFreq,recChannels=recChannels, SNRThresh=SNRThresh, 
                           numUsedFreqs=numUsedFreqs, numPeriodAverages=numPeriodAverages, 
                           numPeriodGrouping=numPeriodGrouping)
  freq = sortFrequencies(freq, bSF, numPeriodGrouping = 1, sortBySNR = sortBySNR)

  if varMeanThresh > 0
    bEmptyTmp = (emptyMeas == nothing) ? bMeas : emptyMeas

    freqVarMean = filterFrequenciesVarMean(bMeas, bEmptyTmp, fgFrames, bgFrames;
                                        thresh=varMeanThresh, minAmplification=minAmplification,
                                        minFreq=minFreq, maxFreq=maxFreq,recChannels=recChannels)
    freq = intersect(freq, freqVarMean)
  end

  if rxNumSamplingPoints(bSF) == rxNumSamplingPoints(bMeas)
    # Ensure that no frequencies are used that are not present in the measurement
    freq = intersect(freq, filterFrequencies(bMeas, numPeriodAverages=numPeriodAverages, 
                                           numPeriodGrouping=numPeriodGrouping))
  end

  @debug "selecting $(length(freq)) frequencies"

  return reconstruction(bSF, bMeas, freq; emptyMeas=emptyMeas, bgFrames=bgFrames, fgFrames=fgFrames, 
                         numPeriodAverages=numPeriodAverages, numPeriodGrouping=numPeriodGrouping, kargs...)
end


function reconstruction(bSF::Union{T,Vector{T}}, bMeas::MPIFile, freq::Array;
  bEmpty = nothing, emptyMeas = bEmpty, bgFrames = 1,
  denoiseWeight = 0, redFactor = 0.0, thresh = 0.0,
  loadasreal = false, solver = Kaczmarz, sparseTrafo = nothing,
  gridsize = gridSizeCommon(bSF), fov=calibFov(bSF), center=[0.0,0.0,0.0], useDFFoV=false,
  deadPixels=nothing, bgCorrectionInternal=false, bgDictSize=nothing, bgFramesDict=nothing,
  numPeriodAverages=1, numPeriodGrouping=1, reco=:default, kargs...) where {T<:MPIFile}

  (typeof(bgFrames) <: AbstractRange && emptyMeas==nothing) && (emptyMeas = bMeas)
  bgCorrection = emptyMeas != nothing ? true : bgCorrectionInternal

  consistenceCheck(bSF, bMeas)

  @debug "Loading System matrix"
  S, grid = getSF(bSF, freq, sparseTrafo, solver;
            bgCorrection=bgCorrection, loadasreal=loadasreal,
            thresh=thresh, redFactor=redFactor,
            useDFFoV=useDFFoV, gridsize=gridsize, fov=fov, center=center,
            deadPixels=deadPixels,numPeriodAverages=numPeriodAverages, 
            numPeriodGrouping=numPeriodGrouping)       

  if denoiseWeight > 0 && sparseTrafo == nothing
    denoiseSF!(S, shape, weight=denoiseWeight)
  end

  # If S is processed and fits not to the measurements because of numPeriodsGrouping
  # or numPeriodAverages being applied we need to set these so that the 
  # measurements are loaded correctly
  if rxNumSamplingPoints(bSF) > rxNumSamplingPoints(bMeas)
    numPeriodGrouping = rxNumSamplingPoints(bSF)  ÷ rxNumSamplingPoints(bMeas)
  end
  if acqNumPeriodsPerFrame(bSF) < acqNumPeriodsPerFrame(bMeas)
    numPeriodAverages = acqNumPeriodsPerFrame(bMeas) ÷ (acqNumPeriodsPerFrame(bSF) * numPeriodGrouping)
  end

  bgDict = getBackgroundDictionary(bSF, bMeas, freq, bgDictSize, bgFramesDict)

  if reco == :default
    return reconstruction(S, bSF, bMeas, freq, grid, emptyMeas=emptyMeas, bgFrames=bgFrames,
                        sparseTrafo=sparseTrafo, loadasreal=loadasreal,
                        bgDict = bgDict,
                        solver=solver, bgCorrectionInternal=bgCorrectionInternal,
                        numPeriodAverages=numPeriodAverages, numPeriodGrouping=numPeriodGrouping; kargs...)
  elseif reco == :tempReg
    return reconstructionTempReg(S, bSF, bMeas, freq, grid, emptyMeas=emptyMeas, bgFrames=bgFrames,
                        sparseTrafo=sparseTrafo, loadasreal=loadasreal,
                        bgDict = bgDict,
                        solver=solver, bgCorrectionInternal=bgCorrectionInternal,
                        numPeriodAverages=numPeriodAverages, numPeriodGrouping=numPeriodGrouping; kargs...)
  else
    error("Parameter :reco is chosen as $(reco), which is not available!")
  end
end

function reconstruction(S, bSF::Union{T,Vector{T}}, bMeas::MPIFile, freq::Array, grid;
  frames = nothing, bEmpty = nothing, emptyMeas= bEmpty, bgFrames = 1, nAverages = 1,
  numAverages=nAverages, bgDict = nothing, bgFramesPost = nothing,
  sparseTrafo = nothing, loadasreal = false, maxload = 100, maskDFFOV=false,
  #weightType=WeightingType.None, weightingLimit = 0, 
  solver = Kaczmarz,
  spectralCleaning=true, spectralLeakageCorrection=spectralCleaning,
  fgFrames=1:10, bgCorrectionInternal=false,
  noiseFreqThresh=0.0, channelWeights=ones(3), 
  numPeriodAverages=1, numPeriodGrouping=1, kargs...) where {T<:MPIFile}

  #(typeof(bgFrames) <: AbstractRange && bEmpty==nothing) && (bEmpty = bMeas)
  bgCorrection = emptyMeas != nothing ? true : false

  @debug "Loading emptymeas ..."
  if emptyMeas!=nothing
    #if acqNumBGFrames(emptyMeas) > 0
    #  uEmpty = getMeasurementsFD(emptyMeas, false, frequencies=freq, frames=bgFrames, #frames=measBGFrameIdx(bEmpty),
    #            numAverages = =length(bgFrames), bgCorrection=bgCorrectionInternal,
    #           loadasreal=loadasreal, spectralLeakageCorrection=spectralLeakageCorrection)
    #else
      uEmpty = getMeasurementsFD(emptyMeas, frequencies=freq, frames=bgFrames, numAverages=length(bgFrames),
      loadasreal = loadasreal,spectralLeakageCorrection=spectralLeakageCorrection, bgCorrection=bgCorrectionInternal,
      numPeriodAverages=numPeriodAverages, numPeriodGrouping=numPeriodGrouping)
      if bgFramesPost != nothing
        uEmptyPost = getMeasurementsFD(emptyMeas, false, frequencies=freq, frames=bgFramesPost,
                    numAverages = length(bgFramesPost), bgCorrection=bgCorrectionInternal,
                    loadasreal = loadasreal, spectralLeakageCorrection=spectralLeakageCorrection,
                    numPeriodAverages=numPeriodAverages, numPeriodGrouping=numPeriodGrouping)
      end
    #end
  end

  frames == nothing && (frames = 1:acqNumFrames(bMeas))

  #=weights = getWeights(weightType, freq, S, weightingLimit=weightingLimit,
                       emptyMeas = emptyMeas, bMeas = bMeas, bgFrames=bgFrames, bSF=bSF,
                       channelWeights = channelWeights)=#

  L = -fld(-length(frames),numAverages) # number of tomograms to be reconstructed
  p = Progress(L, dt=1, desc="Reconstructing data...")

  # initialize sparseTrafo
  B = isnothing(sparseTrafo) ? nothing : createLinearOperator(sparseTrafo, eltype(S), shape = Tuple(shape(grid)))
  @debug "S: $(eltype(S))"
  @debug "B: $(eltype(B))"

  #initialize output
  image = initImage(bSF,bMeas,L,numAverages,grid,false)

  currentIndex = 1
  iterator = numAverages == 1 ? Iterators.partition(frames,maxload) : Iterators.partition(frames,numAverages*maxload)
  for partframes in iterator
    @debug "Loading measurements ..."
    u = getMeasurementsFD(bMeas, frequencies=freq, frames=partframes, numAverages=numAverages,
                          loadasreal=loadasreal, spectralLeakageCorrection=spectralLeakageCorrection,
                          bgCorrection=bgCorrectionInternal, numPeriodAverages=numPeriodAverages, numPeriodGrouping=numPeriodGrouping)
    if emptyMeas!=nothing
      if bgFramesPost == nothing
        u = u .- uEmpty
      else
        for l=1:length(partframes)
          alpha = (partframes[l] - mean(bgFrames)) / (mean(bgFramesPost) - mean(bgFrames))
          u[:,:,l] .-=  (1-alpha).*uEmpty[:,:,1] .+ alpha.*uEmptyPost[:,:,1]
        end
      end
    end

    noiseFreqThresh > 0 && setNoiseFreqToZero(u, freq, noiseFreqThresh, bEmpty = emptyMeas, bMeas = bMeas, bgFrames=bgFrames)

    # convert measurement data if neccessary
    if eltype(S)!=eltype(u)
      @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
      u = map(eltype(S),u)
    end
    @debug "Reconstruction ..."
    c = reconstruction(S, u, bgDict; sparseTrafo=B, progress=p,
		       #weights=weights, 
           solver=solver, shape=shape(grid), kargs...)

    currentIndex = writePartToImage!(image, c, currentIndex, partframes, numAverages)
  end

  return image
end

function writePartToImage!(image, c, currentIndex::Int, partframes, numAverages)
  # permute c's dimensions into image order
  colorsize = size(image,1)
  spatialsize = size(image,2)*size(image,3)*size(image,4)
  inc = -fld(-length(partframes),numAverages)
  c = reshape(c,spatialsize,colorsize,inc)
  c = permutedims(c,[2,1,3])
  # write c to image
  image[Axis{:time}(currentIndex:currentIndex+inc-1)] = c[:]
  currentIndex += inc
  return currentIndex
end

function initImage(bSFs::Union{T,Vector{T}}, bMeas::S, L::Int, numAverages::Int,
		   grid::RegularGridPositions, loadOnlineParams=false) where {T,S<:MPIFile}

  # the number of channels is determined by the number of system matrices
  if isa(bSFs,AbstractVector) || isa(bSFs,MultiContrastFile)
    numcolors = length(bSFs)
    bSF = bSFs[1]
  else
    numcolors = 1
    bSF = bSFs
  end
  # calculate axis
  shp = shape(grid)
  pixspacing, offset = calcSpacingAndOffset(bSF, bMeas, grid)
  dtframes = acqNumAverages(bMeas)*dfCycle(bMeas)*numAverages*1u"s"
  # initialize raw array
  array = Array{Float32}(undef, numcolors,shp...,L)
  # create image
  im = makeAxisArray(array, pixspacing, offset, dtframes)
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
function reconstruction(S, u::Array, bgDict::Nothing=nothing; sparseTrafo = nothing,
                        lambd=0.0, lambda=lambd, λ=lambda, progress=nothing, solver::Type{<:AbstractLinearSolver} = Kaczmarz,
                        weights=nothing, enforceReal=true, enforcePositive=true,
                        relativeLambda=true, reg::Union{Vector{<:AbstractRegularization}, Nothing} = nothing, kargs...)
  N = size(S,2) #prod(shape)
  M = div(length(S), N)

  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N,L)

  norm = NoNormalization()
  if (!isnothing(reg) || sum(abs.(λ)) > 0) && relativeLambda
    norm = SystemMatrixBasedNormalization()
  end
  solverType = eval(Symbol(solver)) # this probably should happen much earlier
  if isnothing(reg)
    reg = AbstractRegularization[L2Regularization(λ)]
    if enforcePositive && !enforceReal
      @warn "enforcePositive also needs enforceReal. Overwriting setting for enforceReal!"
      enforceReal = true
    end
    if enforceReal
      append!(reg, RealRegularization())
    end
    if enforcePositive
      append!(reg, PositiveRegularization())
    end
  else
    if sum(abs.(λ)) > 0
      error("Only λ or an explicit regularization can be given at the same time!")
    end
    if ((RealRegularization() in reg) != enforceReal) || ((PositiveRegularization() in reg) != enforcePositive)
      @warn "An explicit regularization has been given, overriding the behaviour of enforcePositive and enforceReal"
    end
  end

  if !isnothing(sparseTrafo)
    reg = map(r -> TransformedRegularization(r, sparseTrafo), reg)
  end
  
  solv = createLinearSolver(solverType, S; weights=weights,
                            sparseTrafo=sparseTrafo, normalizeReg = norm, reg = reg, kargs...)
  progress==nothing ? p = Progress(L, 1, "Reconstructing data...") : p = progress
  for l=1:L

    d = solve!(solv, u[:,l])

    if !isnothing(sparseTrafo)
      d[:] = sparseTrafo*d #backtrafo from dual space
    end

    #if typeof(B)==LinearSolver.DSTOperator
    #	d=onGridReverse(d,shape)
    #end
    c[:,l] = real( d ) # this one is allocating
    next!(p)
  end

  return c
end

# old code
#if reshapesolution
#  c = reshape(c, shape..., L)
#end
#shape(grid)
