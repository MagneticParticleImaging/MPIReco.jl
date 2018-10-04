using MPIReco
using Test
using Statistics
using Winston

bSF = MPIFile("SF_MP")
b = MPIFile("dataMP01")

c1 = reconstruction(bSF, b;
                   SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2,iterations=1, spectralLeakageCorrection=true)

p = imagesc(c1.data.data[1,:,:,1,1])
savefig(p, "./img/Regular1.png")

# fused lasso
c2 = reconstruction(bSF, b; SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2,iterations=1, solver="fusedlasso", loadasreal=true, lambdaTV=0.1, lambdaL1=0.1)

p = imagesc(c2.data.data[1,:,:,1,1])
savefig(p, "./img/Regular2.png")

# with interpolation
c3 = reconstruction(bSF, b; SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2, gridsize=[100,100,1],iterations=1)

p = imagesc(c3.data.data[1,:,:,1,1])
savefig(p, "./img/Regular3.png")

# with fov adpation and center shift
c4 = reconstruction(bSF, b; SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2, gridsize=[50,50,1],fov=[0.04,0.04,1],
                   center=[0.0,-0.01,0], iterations=1)

p = imagesc(c4.data.data[1,:,:,1,1])
savefig(p, "./img/Regular4.png")


# multi colored reconstruction
c5 = reconstruction([bSF,bSF], b;
                  SNRThresh=5, frames=1, minFreq=80e3,
                  recChannels=1:2,iterations=1)

p = imagesc(c5.data.data[1,:,:,1,1])
savefig(p, "./img/Regular5.png")

p = imagesc(c5.data.data[2,:,:,1,1])
savefig(p, "./img/Regular6.png")


# dict based reco
r = defaultRecoParams()
r[:measPath] = filepath(b)
r[:SFPath] = filepath(bSF)
r[:frames] =  1
r[:minFreq] = 80e3
r[:SNRThresh] = 4
r[:lambd] = 0.001
r[:iterations] = 1

c6 = reconstruction(r)

p = imagesc(c6.data.data[1,:,:,1,1])
savefig(p, "./img/Regular7.png")


###########  reconstruction without storage  ###########

# this is an "in memor
