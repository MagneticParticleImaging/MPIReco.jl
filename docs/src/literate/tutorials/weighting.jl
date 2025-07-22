# # Weighting
include("../../download.jl") #hide

# Often time it is benefical to consider a [weighted least squares problem](https://iopscience.iop.org/article/10.1088/0031-9155/55/6/003) of the form:
# ```math
# \begin{equation}
#   \underset{\mathbf{x}}{argmin} \frac{1}{2}\vert\vert \mathbf{S}\mathbf{c}-\mathbf{u} \vert\vert^2_{\mathbf{W}} + \mathbf{R(u)} .
# \end{equation}
# ```
# where $\mathbf{W}$ is a symmetric, positive weighting matrix and $\vert\vert\mathbf{y}\vert\vert^2_\mathbf{W}$ denotes the weighted Euclidean norm.

# MPIReco provides several different weighting strategies and new strategies can easily be added and plugged into existing algorithms.
# An example of how to implement such a strategy is shown in the How-Tos.
# For an overview of the available strategy consult the API reference.

using MPIReco #hide
bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf")) #hide
b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf")) #hide
params = Dict{Symbol, Any}()
params[:SNRThresh] = 5
params[:frames] = 1:1
params[:minFreq] = 80e3
params[:recChannels] = 1:2
params[:spectralLeakageCorrection] = true
params[:sf] = bSF
params[:reg] = [L2Regularization(0.1f0)];

# To apply different strategies, we simply swap out the `weightingParams` of our single-patch algorithm.
cChannel = reconstruct("SinglePatch", b; params..., weightingParams = ChannelWeightingParameters([0.8, 0.2]))
cRow = reconstruct("SinglePatch", b; params..., weightingParams = RowNormWeightingParameters())
cWhite = reconstruct("SinglePatch", b; params..., weightingParams = WhiteningWeightingParameters(whiteningMeas = bSF));

# We will again visualize the reconstructions with CairoMakie:
using CairoMakie #hide
fig = Figure();
hidedecorations!(heatmap(fig[1, 1], cChannel[1, :, :, 1, 1].data.data, axis = (title = "Channel",)).axis)
hidedecorations!(heatmap(fig[1, 2], cRow[1, :, :, 1, 1].data.data, axis = (title = "Row Norm",)).axis)
hidedecorations!(heatmap(fig[1, 3], cWhite[1, :, :, 1, 1].data.data, axis = (title = "Whitening",)).axis)
fig