using Pkg

# Install required packages
for P in ["HTTP", "PyPlot"]
  !haskey(Pkg.installed(), P) && Pkg.add(P)
end

useCompressedMatrices = true
suffixSM = useCompressedMatrices ? "Small" : "Large"
useFastData = true
suffixMeas = useFastData ? "Fast" : "Slow"

# Download data
include("downloadData.jl")

using MPIReco, PyPlot

################
## Parameters ##
################
Δt = 1.2

# Choose the peak number according to the expected approximate value. Be careful, if there are peaks not corresponding to object motion, division by choosePeak can be false.
# For experiments with different periodic motions it is used to pick one motion frequency (then don't divide by choosePeak)

choosePeak = 4
windowType = 1  # 1: Hann, 2: FT1A05, 3: Rectangle

# Measurement data
datadirMeas = "./data/"
bMeas = MPIFile(datadirMeas*"meas$(suffixMeas).mdf") # high frequency
bBG = MPIFile(datadirMeas*"measBG.mdf") # background measurement

# System matrices
datadirSF = "./data/"
SFall = ["SF$(d)$(suffixSM).mdf" for d=1:4]
bSF = MultiMPIFile(datadirSF.*SFall)

# Background frames
frBG = [1:200, 201:400, 401:600, 601:800]

# Selected frames
recoFrame = 1

####################
## Reconstruction ##
####################

# Reconstruction parameters
lambda = 0.01
iterations = 1

c = reconstructionPeriodicMotion(bSF, bMeas,
                              bEmpty = bBG, frBG = frBG, #background measurement
                              Δt = Δt, choosePeak = choosePeak, recoFrame = recoFrame,
                              windowType=windowType, higherHarmonic=choosePeak,
                              lambda=lambda, iterations=iterations, # reconstruction parameter
                              SNRThresh=2,minFreq=80e3,recChannels=[1,2,3]) #frequency selection

figure(1)
imshow( maximum(arraydata(data(c[1,:,:,:,1])),dims=3)[:,:,1] )
