# # Low-Level Reconstruction
include("../../download.jl") #hide
using MPIReco #hide
bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf")) #hide
b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf")) #hide

# In low-level reconstruction we manually supply both the system matrix $\mathbf{S}$ and the preprocessed measurement vector $\mathbf{u}$.
# As an example we will reproduce a high-level reconstruction with the low-level interface:
params = Dict{Symbol, Any}()
params[:SNRThresh] = 5
params[:sf] = bSF
params[:frames] = 1:acqNumFrames(b)
params[:minFreq] = 80e3
params[:recChannels] = 1:rxNumChannels(b)
params[:iterations] = 1
params[:spectralLeakageCorrection] = true;
params[:reg] = [L2Regularization(0.0f0)]


cHigh = reconstruct("SinglePatch", b; params...);
# While it is possible to manually invoke the processing steps involved in an algorithms reconstruction, we will instead call the respective procsesing functions of
# MPIFiles to prepare $\mathbf{u}$:
freqs = filterFrequencies(bSF, SNRThresh = 5, minFreq = 80e3, recChannels = 1:rxNumChannels(b))
u = getMeasurementsFD(b, frames = 1:acqNumFrames(b), frequencies = freqs, spectralLeakageCorrection=true)

# And afterwards use MPIRecos utility functions to prepare $\mathbf{S}$:
sparseTrafo = nothing
S, grid = getSF(bSF, freqs, sparseTrafo, Kaczmarz)
typeof.([S, grid])

# Now we can configur a low-level reconstruction:
cLow = reconstruct("LowLevel", u; S = S, iterations = params[:iterations], reg = params[:reg])
# Note that the low-level reconstruction returns a matrix without any metadata unlike the other reconstructions.
# The second dimension of the result matrix are the frames. To compare and plot our data we have to reshape it:
sliceLow = sliceLow = reshape(cLow[:, 1], Tuple(grid.shape))
using CairoMakie #hide
fig = Figure();
hidedecorations!(heatmap(fig[1, 1], cHigh[1, :, :, 1, 1].data.data, axis = (title = "High",)).axis)
hidedecorations!(heatmap(fig[1, 2], sliceLow[:, :, 1], axis = (title = "Low",)).axis)
fig
