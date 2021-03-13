using Pkg

# Install required packages
for P in ["HTTP", "PyPlot"]
  !haskey(Pkg.installed(), P) && Pkg.add(P)
end

# Download data
include("downloadData.jl")

using PyPlot, MPIReco

################
## Parameters ##
################

# Measurement data
datadirMeas = "./data/"
f = MPIFile(datadirMeas*"meas.mdf")

# System matrices
datadirSF = "./data/"
fSF = MPIFile(datadirSF*"SF.mdf")

# Background frames
bgFrames = 51:55

# Selected frames
frames = 131  # 51:191

# Bounds for frequency selection
minFreq = 20e3 # change this to 40e3 to exclude the excitation frequency
maxFreq = 200e3

# Use only the first receive channel
recChannels = 1:1

# Solver parameters
λ = 1.0
iterations = 20

# Parameters for the joint estimation of BG and particle concentration
beta = 1.0*5.0.^(-8)
dictSize = 10

# Put all parameters into a dictionary
params = Dict{Symbol,Any}(
  :frames => frames,
  :bgFrames => bgFrames,
  :bgCorrectionInternal => true,
  :λ => λ,
  :recChannels => recChannels,
  :iterations => iterations,
  :minFreq => minFreq,
  :maxFreq => maxFreq,
  :SNRThresh => -1,
  :enforceReal => true,
  :enforcePositive => true,
  :spectralLeakageCorrection =>false
 )

############################
## Perform reconstruction ##
############################

@info "Reconstruction with static BG correction"
paramsStatic = deepcopy(params)
cStaticCorr = reconstruction(fSF,f; paramsStatic...).data.data

@info "Reconstruction with linear interpolated BG correction"
paramsInterp = deepcopy(params)
paramsInterp[:bgFramesPost] = 201:205
cDynCorr = reconstruction(fSF,f; paramsInterp...).data.data

@info "Reconstruction with joint estimation of BG and particle distribution"
paramsJoint = deepcopy(params)
paramsJoint[:bgDictSize] = dictSize
paramsJoint[:beta] = beta
cBGEst = reconstruction(fSF,f; paramsJoint...).data.data

###########################
## Visualize the results ##
###########################

# Windowing parameters
vmin = 0.0
vmax = 0.0006

figure(2, figsize=(6,2))
subplot(1,3,1)
imshow(cStaticCorr[1,:,:,1,1], vmin=vmin, vmax=vmax)
title("Static")
subplot(1,3,2)
imshow(cDynCorr[1,:,:,1,1], vmin=vmin, vmax=vmax)
title("Linear")
subplot(1,3,3)
imshow(cBGEst[1,:,:,1,1], vmin=vmin, vmax=vmax)
title("Joint")
subplots_adjust(left=0.05,bottom=0.05,wspace=0.3,hspace=0.2,right=0.99,top=0.95)
