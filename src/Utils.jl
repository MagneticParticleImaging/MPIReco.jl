export consistenceCheck
import Base: ndims, squeeze

squeeze(A) = squeeze(A,tuple(find(([size(A)...].==1))...))

function generateHeaderDict(bSF::MPIFile, b::MPIFile)
  #header["spatialorder"] TODO
  header = loadMetadata(b)
  header["datatype"] = "MPI"
  header["dim"] = ndims(b)
  header["size"] = calibSize(bSF)
  return header
end


const onlineParams = [:version, :uuid, :time, :dfStrength, :acqGradient,
  :scannerFacility, :scannerOperator, :scannerManufacturer, :scannerName,
  :scannerTopology, :acqFramePeriod, :acqNumPeriodsPerFrame, :acqNumAverages,
  :acqStartTime, :acqOffsetField, :acqOffsetFieldShift,
  :dfNumChannels, :dfPhase, :dfBaseFrequency, :dfDivider,
  :dfCycle, :dfWaveform, :rxNumChannels, :rxBandwidth,
  :rxNumSamplingPoints, :rxTransferFunction] # needs to be extended
  # :tracerName, :tracerBatch, :tracerVendor, :tracerVolume, :tracerConcentration,
  # :tracerSolute, :tracerInjectionTime,
function generateHeaderDictOnline(bSF::MPIFile, b::MPIFile)
   header = loadMetadata(b, onlineParams)
   header["datatype"] = "MPI"
   header["dim"] = ndims(b)
   header["size"] = calibSize(bSF)
  return header
end

# multi patch fallback
generateHeaderDict(bSF::MPIFile, b::Vector) = generateHeaderDict(bSF,b[1])

function consistenceCheck(bSF::MPIFile, bMeas::MPIFile)
  gSF = acqGradient(bSF)[:,:,1,1]
  gMeas = acqGradient(bMeas)[:,:,1,1]
  if gSF != gMeas
    warn("The gradient strength of the system matrix ($gSF T/m) does not fit to the measurements ($gMeas T/m)!")
  end

  dfSF = dfStrength(bSF)
  dfMeas = dfStrength(bMeas)
  if dfSF[:,1] != dfMeas[:,1]
    warn("The drive-field strength of the system matrix ($dfSF mT) does not fit to the measurements ($dfMeas mT)!")
  end

end

function consistenceCheck{T<:MPIFile}(bSFs::Vector{T}, bMeas::MPIFile)
  for bSF in bSFs
    consistenceCheck(bSF,bMeas)
  end
end

function consistenceCheck{T<:MPIFile}(bSF::MPIFile, bMeass::Vector{T})
  for bMeas in bMeass
    consistenceCheck(bSF,bMeas)
  end
end
function consistenceCheck{T<:MPIFile}(bSFs::Vector{T}, bMeass::Vector{T})
  for i = 1:length(bMeass)
    bMeas=bMeass[i]
    bSF=bSFs[i]
    consistenceCheck(bSF,bMeas)
  end
end

#voxelSize(b::MPIFile) = calibFov(b) ./ calibSize(b)
#voxelVolume(b::MPIFile) = prod( voxelSize(b) ) * 1000 #in Liter
#dfFov(b::MPIFile) = squeeze(acqFov(b))
#ffPos(b::MPIFile; kargs...) = acqOffsetFieldShift(b)[:,1,:]
#ndims(b::MPIFile) = sum( (dfStrength(b) .> 0.00001) )
#gridSizeCommon{T<:MPIFile}(bs::Vector{T}) = calibSize(bs[1])
#gridSizeCommon(bs::MPIFile) = calibSize(bs)

###
acqDate(b) = acqStartTime(b)
subjectName(b) = experimentSubject(b)
scanName(b) = experimentName(b)
expno(b) = experimentNumber(b)
numAverages(b) = acqNumAverages(b)
frequencies(b) = rxFrequencies(b)
operator(b) = scannerOperator(b)
numReceivers(b) = rxNumChannels(b)
bandwidth(b) = rxBandwidth(b)
tracer(b) = tracerName(b)
function description(b)
   println(filepath(b))
   experimentDescription(b)
end
getSNRAllFrequencies(b) = calibSNR(b)[:,:,1]
measPath(b::BrukerFile) = b.path
measPath(b::MDFFile) = b.filename
gridSize(b::MPIFile) = squeeze(calibSize(b))
measPath{T<:MPIFile}(bs::Vector{T}) = [measPath(b) for b in bs]
gridSize{T<:MPIFile}(bs::Vector{T}) = [gridSize(b) for b in bs]
gridSizeCommon{T<:MPIFile}(bs::Vector{T}) = gridSize(bs[1])
gridSizeCommon(bs::MPIFile) = gridSize(bs)
fov(b::MPIFile) = calibFov(b)
fov{T<:MPIFile}(b::Vector{T}) = fov(b[1])
sfGradient(b::MPIFile) = diag(acqGradient(b)[:,:,1,1])
sfGradient(b::MPIFile,dim) = sfGradient(b)[dim]
numFreq(b::MPIFile) = rxNumFrequencies(b)
ffPos(b::MPIFile; kargs...) = squeeze(acqOffsetFieldShift(b))
function ffPos{T<:BrukerFile}(b::Vector{T};alpha=[0,0,0])
  fovCenter = zeros(3,length(b))
  for l=1:length(b)
    fovCenter[:,l] = ffPos(b[l], alpha=alpha)
  end
  return fovCenter
end
numScans(b) = acqNumFrames(b)
dfcycle(b::MPIFile) = dfCycle(b)
ndims(b::MPIFile) = sum( (dfStrength(b) .> 0.00001) )
numPatches(b) = acqNumPeriodsPerFrame(b)
voxelSize(b::MPIFile) = fov(b) ./ gridSize(b)
voxelVolume(b::MPIFile) = prod( voxelSize(b) ) * 1000 #in Liter
dfFov(b) = squeeze(acqFov(b))
numTimePoints(b::MPIFile) = rxNumSamplingPoints(b)
