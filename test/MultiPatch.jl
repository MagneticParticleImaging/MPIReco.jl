using MPIReco
using Base.Test

using PyPlot


# MultiPatch
bSF = MPIFile("SF_MP")

b = MultiMPIFile(["dataMP01", "dataMP02", "dataMP03", "dataMP04"])
figure(7)

@time c1 = reconstruction(bSF, b;
                   SNRThresh=5, frames=1, minFreq=80e3, OffsetSF=0,OffsetFF=[0.008,0.008,0.0005],
                   recChannels=1:2,iterations=1, spectralLeakageCorrection=false)

subplot(2,2,1)
imshow(c1[:,:,1,1])

@time c2 = reconstruction(bSF, b;
                   SNRThresh=5, frames=1, minFreq=80e3, OffsetSF=0,OffsetFF=[0.008,0.008,0.0005],
                   recChannels=1:2,iterations=1, roundPatches=true, spectralLeakageCorrection=false)

subplot(2,2,2)
imshow(c2[:,:,1,1])



# MultiPatch MultiSF
bSFs = MultiMPIFile(["SF_MP01", "SF_MP02", "SF_MP03", "SF_MP04"])

@time c3 = reconstruction(bSFs, b;
                   SNRThresh=5, frames=1, minFreq=80e3, OffsetSF=0,OffsetFF=[0.0001,0.0001,0.0005],
                   recChannels=1:2,iterations=1, spectralLeakageCorrection=false)

subplot(2,2,3)
imshow(c3[:,:,1,1])
