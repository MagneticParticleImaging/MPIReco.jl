



function reconstructionTempReg(S, bSF::Union{T,Vector{T}}, bMeas::MPIFile, freq::Array, grid;
    frames = nothing, bEmpty = nothing, emptyMeas= bEmpty, bgFrames = 1, nAverages = 1,
    numAverages=nAverages, bgDict = nothing, 
    sparseTrafo = nothing, loadasreal = false, maxload = 100, maskDFFOV=false,
    weightType=WeightingType.None, weightingLimit = 0, solver = "kaczmarz",
    spectralCleaning=true, spectralLeakageCorrection=spectralCleaning,
    fgFrames=1:10, bgCorrectionInternal=false,
    noiseFreqThresh=0.0, channelWeights=ones(3), 
    numPeriodAverages=1, numPeriodGrouping=1, kargs...) where {T<:MPIFile}
  
    bgCorrection = emptyMeas != nothing ? true : false
  
    @debug "Loading emptymeas ..."
    if emptyMeas!=nothing
        uEmpty = getMeasurementsFD(emptyMeas, frequencies=freq, frames=bgFrames, numAverages=length(bgFrames),
        loadasreal = loadasreal,spectralLeakageCorrection=spectralLeakageCorrection, bgCorrection=bgCorrectionInternal,
        numPeriodAverages=numPeriodAverages, numPeriodGrouping=numPeriodGrouping)
    end
  
    frames == nothing && (frames = 1:acqNumFrames(bMeas))
  
    weights = getWeights(weightType, freq, S, weightingLimit=weightingLimit,
                         emptyMeas = emptyMeas, bMeas = bMeas, bgFrames=bgFrames, bSF=bSF,
                         channelWeights = channelWeights)
  
    L = -fld(-length(frames),numAverages) # number of tomograms to be reconstructed
    p = Progress(L, 1, "Reconstructing data...")
  
    # initialize sparseTrafo
    B = linearOperator(sparseTrafo, shape(grid), eltype(S))
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
                 weights=weights, solver=solver, shape=shape(grid), kargs...)
  
      currentIndex = writePartToImage!(image, c, currentIndex, partframes, numAverages)
    end
  
    return image
  end



mutable struct TemporalRegularizationOperator{V<:AbstractMatrix, U<:Positions}
    S::V
    Φ::V
    grid::U
    N::Int
    M::Int
    L::Int # Number of all raw data frames
    Θ::Int # Number of subsampled foreground frames 
    Γ::Int # Number of subsampled background frames 
    θ::Vector{Int} # Sampling indices (sorted)
    γ::Vector{Int} # Sampling indices (sorted)
  end
  
  eltype(TempOp::TemporalRegularizationOperator) = eltype(TempOp.S)