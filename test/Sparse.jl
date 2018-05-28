using MPIReco
using Base.Test

using PyPlot



bSF = MPIFile("SF_MP")
b = MPIFile("dataMP01")

c1 = reconstruction(bSF, b; lambd=0.1,
                   SNRThresh=10, frames=1, minFreq=80e3, sparseTrafo="FFT",
                   recChannels=1:2,iterations=10, spectralLeakageCorrection=true,
                   redFactor =0.1)

figure(1)
imshow(c1[1,:,:,1,1])
