using MPIReco
using Base.Test

using PyPlot



bSF = MPIFile("SF_MP")

b = MPIFile("dataMP01")

c = reconstruction(bSF, b;
                   SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2,iterations=1, spectralLeakageCorrection=true)


figure(1)
imshow(c[1,:,:,1,1])
