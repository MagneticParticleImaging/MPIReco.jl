using MPIReco
using Test
using Statistics
using Winston


# MultiPatch
bSF = MultiMPIFile(["SF_MP"])

b = MultiMPIFile(["dataMP01", "dataMP02", "dataMP03", "dataMP04"])

@time c1 = reconstruction(bSF, b;
                   SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2,iterations=1, spectralLeakageCorrection=false)

#figure(11)
imagesc(c1[1,:,:,1,1])

@time c2 = reconstruction(bSF, b;
                   SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2,iterations=1, roundPatches=true, spectralLeakageCorrection=false)

#figure(12)
imagesc(c2[1,:,:,1,1])



# MultiPatch MultiSF
bSFs = MultiMPIFile(["SF_MP01", "SF_MP02", "SF_MP03", "SF_MP04"])

@time c3 = reconstruction(bSFs, b;
                   SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2,iterations=1, spectralLeakageCorrection=false)

#figure(12)
imagesc(c3[1,:,:,1,1])



# Flexible MultiPatch construction
bSFs = MultiMPIFile(["SF_MP01", "SF_MP02", "SF_MP03", "SF_MP04"])
mapping = [1,2,3,4]
freq = filterFrequencies(bSFs, SNRThresh=5, minFreq=80e3)
S = [getSF(SF,freq,nothing,"kaczmarz", bgcorrection=false)[1] for SF in bSFs]
SFGridCenter = zeros(3,4)
FFPos = zeros(3,4)
FFPos[:,1] = [-0.008, 0.008, 0.0]
FFPos[:,2] = [-0.008, -0.008, 0.0]
FFPos[:,3] = [0.008, 0.008, 0.0]
FFPos[:,4] = [0.008, -0.008, 0.0]

@time c4 = reconstruction(bSFs, b;
                   SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2,iterations=1, spectralLeakageCorrection=false, mapping=mapping,
                   systemMatrices = S, SFGridCenter=SFGridCenter, FFPos=FFPos, FFPosSF=FFPos)

#figure(13)
imagesc(c4[1,:,:,1,1])
