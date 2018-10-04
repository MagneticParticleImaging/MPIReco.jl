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

p = imagesc(c1.data.data[1,:,:,1,1])
savefig(p, "./img/Sparse1.png")


c2 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="FFT",
                   recChannels =1:2,
                   iterations = 3,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = false)

p = imagesc(c2.data.data[1,:,:,1,1])
savefig(p, "./img/Sparse2.png")

c3 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="FFT",
                   recChannels =1:2,
                   iterations = 3,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = true)

p = imagesc(c3.data.data[1,:,:,1,1])
savefig(p, "./img/Sparse3.png")

c2 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DCT",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = false)

p = imagesc(c2.data.data[1,:,:,1,1])
savefig(p, "./img/Sparse4.png")

c3 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DCT",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=true,
                   redFactor = redFactor,
                   useDFFoV = true)

p = imagesc(c3.data.data[1,:,:,1,1])
savefig(p, "./img/Sparse5.png")

c2 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DST",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=false,
                   redFactor = redFactor,
                   useDFFoV = false)

p = imagesc(c2.data.data[1,:,:,1,1])
savefig(p, "./img/Sparse6.png")


c3 = reconstruction(bSF, b; lambd=0.01,
                   SNRThresh = 3, frames=1:100, minFreq=80e3, nAverages=100,
                   sparseTrafo="DST",
                   recChannels =1:2,
                   iterations = 1,
                   spectralLeakageCorrection=true,
                   redFactor = redFactor,
                   useDFFoV = true)

p = imagesc(c3.data.data[1,:,:,1,1])
savefig(p, "./img/Sparse7.png")
