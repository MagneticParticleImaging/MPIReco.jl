using MPIReco
using Base.Test

using PyPlot


# MultiPatch
bSF = MPIFile("SF_MP")

b = MultiMPIFile(["dataMP01", "dataMP02", "dataMP03", "dataMP04"])

c = reconstruction(bSF, b;
                   SNRThresh=5, frames=1, minFreq=80e3, OffsetSF=0,OffsetFF=[0.0001,0.0001,0.0005],
                   recChannels=1:2,iterations=1, spectralLeakageCorrection=false)


figure(2)
imshow(c[:,:,1,1])
