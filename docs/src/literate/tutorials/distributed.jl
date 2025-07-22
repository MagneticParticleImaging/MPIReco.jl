# # Distributed Image Reconstruction
include("../../download.jl") #hide
# MPIReco image reconstructions can also be executed across different computers, thanks to [DaggerImageReconstruction](https://github.com/JuliaImageRecon/DaggerImageReconstruction.jl).
# With this feature one could perform a reconstruction over the network on a machine with some specific resource, be it a specific system matrix or access to a GPU.

# To enable a distributed reconstruction one has to first use the Distributed standard library to add a new Julia process on the remote machine.
# ```julia
# using Distributed
# workers = addprocs(["remote_address])
# worker = first(workers)
# ```
using Distributed #hide
worker = 1 #hide # comment out to properly connect to a different process

# The worker is the Julia process id, which is used to identify where to move data to and perform computations on.
# Afterwards we just load both MPIReco and DaggerImageReconstruction:
using MPIReco, DaggerImageReconstruction

# Instead of a MPIFile, we now want to create a DMPIFile, a distributed MPIFIle. Such a file expects a path on the remote machine together with the worker id:
distr_bSF = DMPIFile(joinpath(datadir, "calibrations", "12.mdf"), worker = worker)
distr_b = DMPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"), worker = worker)

# Note that you might need to consider differences between the operating systems of both machines.
# For example, in this case we constructing our filepath locally, while evaluating it on the remote.
# You could also construct the path on the worker with:
remotecall_fetch(() -> joinpath("..."), worker)

# Once we have our distributed files, we can configure our reconstruction like usual and
# the algorithm figures out on which worker to execute based on the provided distributed MPIFiles.
# In this case, both files were located on the remote. If the measurements are located on the local machine, then one has to transfer the data over the network.
c1 = reconstruct("SinglePatch", distr_b;
                   SNRThresh=5,
                   sf = distr_bSF,
                   frames=1:acqNumFrames(distr_b),
                   minFreq=80e3,
                   recChannels=1:rxNumChannels(distr_b),
                   iterations=1,
                   spectralLeakageCorrection=true);

# Most MPIFiles instances are just handles to one or more local files and thus can't be meaningfully send over the network.
# To get around this, you can transform an MDF into an MDFinMemory, which is a fully in-memory representation of an MDF.
b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"))
bInMemory = MDFv2InMemory(b)

# Now the algorithm can't determinte the worker from the measurement anymore and we have to utilize the `loadDaggerPlan` function from DaggerImageReconstruction.
# This function expects a path to a local RecoPlan, which is then constructed on the given worker:
planpaths = getRecoPlanList(;full=true)
index = findfirst(endswith("SinglePatch.toml"), planpaths)
distr_plan = loadDaggerPlan(planpaths[index], getRecoPlanModules(), worker = 1)

# We can then configure the plan as usual:
setAll!(distr_plan; SNRThresh=5,
                   sf = distr_bSF,
                   frames=1:acqNumFrames(bInMemory),
                   minFreq=80e3,
                   recChannels=1:rxNumChannels(bInMemory),
                   iterations=1,
                   spectralLeakageCorrection=true)

# And perform a normal reconstruction with it:
distr_algo = build(distr_plan)
c2 = reconstruct(distr_algo, bInMemory)
isapprox(c1.data.data, c2.data.data)