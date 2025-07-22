# # Getting Started

# To begin, we first need to gather some MPI data. To do this we could use data from the [OpenMPIData](https://github.com/MagneticParticleImaging/OpenMPIData.jl) initiative, however for these examples we will use the data used in testing MPIReco.

# To access this test data, first, enter the Pkg mode in Julia (`]`) and execute the unit tests of MPIReco::
# ```julia
# test MPIReco
# ```
# This will download and unpack some MPI measurements and calibration MDF files and perform tests with the data.
# The download location will be printed at the start of the test execution. You can cancel the test execution 
# once the data is downloaded by pressing `Ctrl + C` or closing the Julia process.

# After the download, several MPI files will be present in the specified directory. 
# All subsequent examples assume that you have assigned the path to this directory to a variable named `datadir`:
include("../../download.jl"); #hide
# ```julia
# const datadir = joinpath("...","artifacts", "...") # enter path here
# ```

# ## First Reconstruction

# We will start looking at a small reconstruction script. First, we load MPIReco:
using MPIReco

# Next, we open handles for the system matrix and measurement data. 
# Both are created via MPIFile function, which can be, for instance, an MDFFile or a BrukerFile.
bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf"))
b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"))

# To interact with the files, you can use the functionality of MPIFiles, which is also exported by MPIReco:
acqNumFrames(b)
# We refer to the documentation for MPIFiles.jl and MDF for more information on the available data and functions.

# Now we will perform a system-matrix based reconstruction using the `reconstruct` function. 
# The function takes as arguments the algorithm to use, the measurement data, and various keyword arguments.
# The available parameters and their interpretation depend on the chosen algorithm. 
# In this case, the algorithm also expects the system matrix to be set. 
# We set the SNR threshold to 5, meaning that only matrix rows with an SNR above 5 will be used during reconstruction.
# The parameter `frame` decides which frame of the measured data should be reconstructed.

c = reconstruct("SinglePatch", b;
                   SNRThresh=5,
                   sf = bSF,
                   frames=1:acqNumFrames(b),
                   minFreq=80e3,
                   recChannels=1:rxNumChannels(b),
                   iterations=1,
                   spectralLeakageCorrection=true)

# Notice that the result is not just an array of particle concentration, but contains a variety of metadata as well.
# This is because `c` is a:
typeof(c)

# You can access metadata as properties, such as:
c.tracerName

# The result can also be treated like a normal 5-dimensional array.
# The first dimension is for different contrasts, followed by the three spatial dimensions 
# and lastly the temporal dimension of the different frames.
size(c)

# You can access the underlying data just like any other array in Julia. 
# For example, you can access an XY slice as follows:
slice = c[1, :, :, 1, 1]
typeof(slice)

# Note that the slice is still an ImageMeta array.
# Inside the array is an AxisArray, another array type with metadata.
# Not every Julia package has methods for such metadata arrays. 
# To visualize our results with CairoMakie, we need to extract the underlying data.
sliceRaw = slice.data.data
typeof(sliceRaw)

# We can now visualize our first MPI reconstruction using MPIReco:
using CairoMakie
fig = heatmap(sliceRaw)
hidedecorations!(fig.axis)
fig