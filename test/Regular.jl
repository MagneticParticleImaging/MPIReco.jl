using MPIReco
using Test
using Statistics
using Winston



bSF = MPIFile("SF_MP")
b = MPIFile("dataMP01")

c1 = reconstruction(bSF, b;
                   SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2,iterations=1, spectralLeakageCorrection=true)

#figure(1)
imagesc(c1[1,:,:,1,1])


# fused lasso
c2 = reconstruction(bSF, b; SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2,iterations=1, solver="fusedlasso", loadasreal=true, lambdaTV=0.1, lambdaL1=0.1)

#figure(2)
imagesc(c2[1,:,:,1,1])

# with interpolation
c3 = reconstruction(bSF, b; SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2, gridsize=[100,100,1],iterations=1)

#figure(3)
imagesc(c3[1,:,:,1,1])

# with fov adpation and center shift
c4 = reconstruction(bSF, b; SNRThresh=5, frames=1, minFreq=80e3,
                   recChannels=1:2, gridsize=[50,50,1],fov=[0.04,0.04,1],
                   center=[0.0,-0.01,0], iterations=1)

#figure(4)
imagesc(c4[1,:,:,1,1])


# multi colored reconstruction
c5 = reconstruction([bSF,bSF], b;
                  SNRThresh=5, frames=1, minFreq=80e3,
                  recChannels=1:2,iterations=1)

#figure(5)
imagesc(c5[1,:,:,1,1])
#figure(6)
imagesc(c5[2,:,:,1,1])


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

#figure(7)
imagesc(c6[1,:,:,1,1])


###########  reconstruction without storage  ###########

# this is an "in memor
