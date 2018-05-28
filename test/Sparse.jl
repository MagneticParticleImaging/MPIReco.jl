using MPIReco
using Base.Test

using PyPlot


bSF = MPIFile("systemMatrix.mdf")
b = MPIFile("measurement.mdf")


c1 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
                   recChannels =1:2,
                   iterations = 3,
                   spectralLeakageCorrection=false)


figure(1)
subplot(2,2,1)
imshow(c1[1,:,:,1,1])


c2 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DCT",
                   recChannels =1:2,
                   iterations = 3,
                   spectralLeakageCorrection=false,
                   redFactor = 0.01,
                   useDFFoV = false)

subplot(2,2,2)
imshow(c2[1,:,:,1,1])


c3 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DCT",
                   recChannels =1:2,
                   iterations = 3,
                   spectralLeakageCorrection=false,
                   redFactor = 0.01,
                   useDFFoV = true)

subplot(2,2,3)
imshow(c3[1,:,:,1,1])
