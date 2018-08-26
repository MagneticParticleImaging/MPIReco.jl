using MPIReco
using Test
using Statistics
using Winston


# MultiGradient

# low gradient with patch no 3
b = MultiMPIFile(["dataMG_G1", "dataMG_G2_03"])
bSFs = MultiMPIFile(["SF_MG_G1", "SF_MG_G2"])

@time c1 = reconstruction(bSFs, b;
                   SNRThresh=2, frames=1, lambd=0.003, minFreq=80e3,
                   recChannels=1:2,iterations=3, roundPatches=false)

figure(8)
imagesc(reverse(c1[1,:,:,1,1]',dims=1))

# low gradient with all 3 high gradient patches
b = MultiMPIFile(["dataMG_G1", "dataMG_G2_01", "dataMG_G2_02", "dataMG_G2_03"])

@time c2 = reconstruction(bSFs, b;
                   SNRThresh=2, frames=1, lambd=0.003, minFreq=80e3,
                   recChannels=1:2,iterations=3, roundPatches=false)

figure(9)
imagesc(reverse(c2[1,:,:,1,1]',dims=1))
