# # Multi-Patch Reconstruction
include("../../download.jl") #hide
# For multi-patch reconstruction the method proposed by [Szwargulski et al.](https://www.ncbi.nlm.nih.gov/pubmed/30334751) is implemented in MPIReco. It is generalized however as described in [Boberg et al.](https://doi.org/10.1109/TMI.2019.2949171).

# We first discuss the measurements for the multi-patch case. On modern MPI
# scanners the `BrukerFile` or `MDFFile` can be used as is. However, the data
# that we use in our unit tests consists of several single-patch measurements.
# We therefore combine them manually into an `MultiMPIFile`:
using MPIReco #hide
files = ["1.mdf", "2.mdf", "3.mdf", "4.mdf"]
b = MultiMPIFile(joinpath.(datadir, "measurements", "20211226_203916_MultiPatch", files))

# And now `b` can be used as if it was a multi-patch file.

# Next we take a look at the system matrix. The most simple approach is to use a single system
# matrix that was measured at the center. This can be done using
bSF = MultiMPIFile([joinpath(datadir, "calibrations", "12.mdf")])

# Afterwards we can perform a reconstruction with 
c1 = reconstruct("MultiPatch", b; sf = bSF, SNRThresh=5, frames=1:acqNumFrames(b), minFreq=80e3,
                   recChannels=1:rxNumChannels(b), iterations=1, spectralLeakageCorrection=false, tfCorrection = false);
# The parameters are essentially the same as in the previous reconstructions, we only change the algorithm and the MPIFiles.

# It is also possible to use multiple system matrices, which is currently the
# best way to take field imperfection into account. Our test data has four patches
# and we therefore can use
sf_files =  ["8.mdf", "9.mdf", "10.mdf", "11.mdf"]
bSFs = MultiMPIFile(joinpath.(datadir, "calibrations", sf_files))

c2 = reconstruct("MultiPatch", b; sf = bSFs, SNRThresh=5, frames=1:acqNumFrames(b), minFreq=80e3,
                   recChannels=1:rxNumChannels(b), iterations=1, spectralLeakageCorrection=false, tfCorrection = false);

# Now we want somewhat more flexibility and
# * define a mapping between the system matrix and the patches, here we allow to
#   use the same system matrix for multiple patches
# * make it possible to change the FFP position. Usually the value stored in the
#   file is not 100% correct due to field imperfections.
# To achieve this, we swap out parts of the MultiPatch blueprint with a parameter which allows us to explicitly set all these parameters

# We perform our own frequency filtering.
freq = filterFrequencies(bSFs, SNRThresh=5, minFreq=80e3);

# And load four different system matrices for each patch.
S = [getSF(SF,freq,nothing,"Kaczmarz", bgcorrection=false, tfCorrection = false)[1] for SF in bSFs]

# Afterwards we can describe our patch mapping and positions:
mapping = [1,2,3,4]
SFGridCenter = zeros(3,4)
FFPos = zeros(3,4)
FFPos[:,1] = [-0.008, 0.008, 0.0]
FFPos[:,2] = [-0.008, -0.008, 0.0]
FFPos[:,3] = [0.008, 0.008, 0.0]
FFPos[:,4] = [0.008, -0.008, 0.0];

# Lastly, we wrap everything in a parameter structs:
opParams = ExplicitMultiPatchParameter(;tfCorrection = false, systemMatrices = S, SFGridCenter = SFGridCenter, mapping = mapping)
ffPos = CustomFocusFieldPositions(FFPos);

# and perform our reconstruction:
c3 = reconstruct("MultiPatch", b; sf = bSFs, opParams = opParams, ffPos = ffPos, ffPosSF = ffPos,
                  SNRThresh=5, frames=1:acqNumFrames(b), minFreq=80e3,
                  recChannels=1:rxNumChannels(b), iterations=1, spectralLeakageCorrection=false, tfCorrection = false);

# We can again visualize our different multi-patch reconstructions:
using CairoMakie #hide
fig = Figure();
hidedecorations!(heatmap(fig[1, 1], c1[1, :, :, 1, 1].data.data, axis = (title = "Single SM",)).axis)
hidedecorations!(heatmap(fig[1, 2], c2[1, :, :, 1, 1].data.data, axis = (title = "Multiple SM",)).axis)
hidedecorations!(heatmap(fig[1, 3], c3[1, :, :, 1, 1].data.data, axis = (title = "Explicit SM",)).axis)
fig