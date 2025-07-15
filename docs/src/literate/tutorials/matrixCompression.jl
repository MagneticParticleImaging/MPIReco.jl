# # Matrix-Compression Techniques
include("../../download.jl") #hide
# Matrix compression can significantly accelerate the reconstruction process.
# This is achieved by transforming the system matrix `S` into a different domain 
# through a basis transformation applied to its rows.

# In `MPIReco.jl`, matrix compression can be enabled by specifying a sparse 
# system-matrix loading parameter and selecting the desired sparsity transformation.
using MPIReco #hide
bSF = MPIFile(joinpath(datadir, "calibrations", "7.mdf")) #hide
b = MPIFile(joinpath(datadir, "measurements","20211226_204612_Dice", "1.mdf")) #hide
params = Dict{Symbol, Any}()
params[:SNRThresh] = 2
params[:frames] = 1:100
params[:numAverages] = 100
params[:minFreq] = 80e3
params[:recChannels] = 1:2
params[:iterations] = 3
params[:spectralLeakageCorrection] = false
params[:sf] = bSF;

# We first reconstruct as usual:
cDense = reconstruct("SinglePatch", b; params...);

# Afterwards we use a blueprint with a sparse system-matrix loading and setting a sparsity transformation:
cSparse = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "FFT");
# Possible transformations are "FFT", "DCT_IV" and "DST".

# The transformations can be restricted to the drive-field field-of-view by setting `useDFFoV = true`. 
# The compression factor, which controls how many coefficients are dropped 
# after applying the transformation, is determined by the `redFactor` parameter.
# For example, a reduction factor of `redFactor = 0.01` will drop 99% of the data.
cRed = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "FFT", useDFFoV = true, redFactor = 0.01);

# We can again visualize the data with CairoMakie:
using CairoMakie #hide
fig = Figure();
hidedecorations!(heatmap(fig[1, 1], cDense[1, :, :, 1, 1].data.data, axis = (title = "Dense",)).axis)
hidedecorations!(heatmap(fig[1, 2], cSparse[1, :, :, 1, 1].data.data, axis = (title = "Sparse",)).axis)
hidedecorations!(heatmap(fig[1, 3], cRed[1, :, :, 1, 1].data.data, axis = (title = "Sparse DF FOV",)).axis)
fig