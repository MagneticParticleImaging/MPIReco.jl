using MPIReco
using Test
using Statistics
using Winston


bSF = MPIFile("systemMatrix.mdf")
b = MPIFile("measurement.mdf")


redFactor = 0.01

c1 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
                   recChannels =1:2,
                   iterations = 3,
                   spectralLeakageCorrection=false)


figure(16)
imagesc(c1[1,:,:,1,1])


c2 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="FFT",
                   recChannels =1:2,
                   iterations = 3,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = false)

figure(17)
imagesc(c2[1,:,:,1,1])


c3 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="FFT",
                   recChannels =1:2,
                   iterations = 3,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = true)

figure(18)
imagesc(c3[1,:,:,1,1])

c2 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DCT",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = false)

figure(19)
imagesc(c2[1,:,:,1,1])


c3 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DCT",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=true,
                   redFactor = redFactor,
                   useDFFoV = true)

figure(20)
imagesc(c3[1,:,:,1,1])

c2 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DST",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = false)

figure(21)
imagesc(c2[1,:,:,1,1])


c3 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DST",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=true,
                   redFactor = redFactor,
                   useDFFoV = true)

figure(22)
imagesc(c3[1,:,:,1,1])
