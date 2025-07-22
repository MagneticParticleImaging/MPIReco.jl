# # Multi-Contrast Reconstruction
include("../../download.jl") #hide

# So far we have discussed single-contrast reconstructions of the form:
# ```math
# \begin{equation}
#   \underset{\mathbf{c}}{argmin} \frac{1}{2}\vert\vert \mathbf{S}\mathbf{c}-\mathbf{u} \vert\vert_2^2 + + \mathbf{R(x)}
# \end{equation}
# ```
# where $\mathbf{S}$ is a single system matrix, $\mathbf{u}$ is the measurement vector, and $\mathbf{R(x)}$ is an (optional) regularization term.
# In a multi-contrast reconstruction one can use two or more system matrices to solve:
# ```math
# \begin{equation}
#    \underset{c1c2}{argmin} \frac{1}{2}\left\| \mathbf{S_1S_2} \begin{pmatrix} \mathbf{c_1} \\ \mathbf{c_2} \end{pmatrix} - \mathbf{u} \right\|_2^2 + \mathbf{R(x)}
# \end{equation}
# ```

# To apply this technique to our algorithms, we simply have to provide multiple system matrices to our reconstructions:
using MPIReco #hide
bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf")) #hide
b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf")) #hide
bContrast = MultiContrastFile([bSF, bSF])

# and can then reconstruct as usual:
c = reconstruct("SinglePatch", b;
                   SNRThresh=5,
                   sf = bContrast,
                   frames=1:acqNumFrames(b),
                   minFreq=80e3,
                   recChannels=1:rxNumChannels(b),
                   iterations=1,
                   spectralLeakageCorrection=true);
size(c)
# Note that `c` now has two entries in its first dimension, one for each contrast. 


# We can then just simply access the data of each contrast by accessing the specific dimensions:
using CairoMakie #hide
fig = Figure();
hidedecorations!(heatmap(fig[1, 1], c[1, :, :, 1, 1].data.data, axis = (title = "Contrast 1",)).axis)
hidedecorations!(heatmap(fig[1, 2], c[2, :, :, 1, 1].data.data, axis = (title = "Contrast 2",)).axis)
fig
# In this case the channels are identical, because we reused the same system matrix.