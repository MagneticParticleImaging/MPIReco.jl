using MPIReco
using Base.Test

using PyPlot


bSF = MPIFile("systemMatrix.mdf")
b = MPIFile("measurement.mdf")


redFactor = 0.01

c1 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
                   recChannels =1:2,
                   iterations = 3,
                   spectralLeakageCorrection=false)


figure(1)
subplot(5,2,1)
imshow(c1[1,:,:,1,1])


c2 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="FFT",
                   recChannels =1:2,
                   iterations = 3,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = false)

subplot(5,2,3)
imshow(c2[1,:,:,1,1])


c3 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="FFT",
                   recChannels =1:2,
                   iterations = 3,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = true)

subplot(5,2,4)
imshow(c3[1,:,:,1,1])

c2 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DCT",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = false)

subplot(5,2,5)
imshow(c2[1,:,:,1,1])


c3 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DCT",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=true,
                   redFactor = redFactor,
                   useDFFoV = true)

subplot(5,2,6)
imshow(c3[1,:,:,1,1])

c2 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DST",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = false)

subplot(5,2,7)
imshow(c2[1,:,:,1,1])


c3 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DST",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=true,
                   redFactor = redFactor,
                   useDFFoV = true)

subplot(5,2,8)
imshow(c3[1,:,:,1,1])
