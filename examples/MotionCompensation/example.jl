using Pkg

# Install required packages
for P in ["HTTP", "PyPlot"]
  !haskey(Pkg.installed(), P) && Pkg.add(P)
end


using HTTP

# Doanload data
if !isdir("data")
  @info "download data.zip"
  HTTP.open("GET", "http://media.tuhh.de/ibi/MotionCompensation/data.zip") do http
    open("data.zip", "w") do file
        write(file, http)
    end
  end
  @info "extracting data.zip"
  run(`unzip -oq data.zip`)
  rm("data.zip")
end

using MPIReco, PyPlot

################
## Parameters ##
################
alpha = 0.4

# Choose the peak number according to the expected approximate value. Be careful, if there are peaks not corresponding to object motion, division by choosePeak can be false.
# For experiments with different periodic motions it is used to pick one motion frequency (then don't divide by choosePeak)

choosePeak = 4

windowType = 1  # 1: Hann, 2: FT1A05, 3: Rectangle

# Measurement data
datadirMeas = "./data/"
bMeas = MPIFile(datadirMeas*"measFast.mdf") # high frequency
bBG = MPIFile(datadirMeas*"measBG.mdf") # background measurement

# System matrices
datadirSF = "./data/"
SFall = ["SF1Small.mdf","SF2Small.mdf","SF3Small.mdf","SF4Small.mdf"]
bSF = MultiMPIFile(datadirSF.*SFall)

# Background frames
frBG = [1:200,201:400, 401:600 ,601:800]

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
                              alpha = alpha, choosePeak = choosePeak, recoFrame = recoFrame,
                              windowType=windowType, higherHarmonic=choosePeak,
                              lambda=lambda, iterations=iterations, # reconstruction parameter
                              SNRThresh=10,minFreq=80e3,recChannels=[1,2,3]) #frequency selection

figure(1)
imshow( maximum(arraydata(data(c[1,:,:,:,1])),dims=3)[:,:,1] )
