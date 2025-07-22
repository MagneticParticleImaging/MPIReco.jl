# # Basic Reconstruction
include("../../download.jl") #hide
# MPIReco.jl provides different reconstruction interfaces. 
# All of the interfaces interact with `RecoPlans` from `AbstractImageReconstruction`.
# These are blueprints from which different algorithms can be constructed. 
# from which different algorithms can be constructed. For more details on these blueprints, we refer to the [AbstractImageReconstruction documentation](https://juliaimagerecon.github.io/AbstractImageReconstruction.jl/dev/).

# The reconstruction algorithm in MPIReco mainly focus on system-matrix based reconstructions, where one considers an inverse problem of the form:
# ```math
# \begin{equation}
#   \underset{\mathbf{c}}{argmin} \frac{1}{2}\vert\vert \mathbf{S}\mathbf{c}-\mathbf{u} \vert\vert_2^2 + + \mathbf{R(x)}
# \end{equation}
# ```
# where $\mathbf{S}$ is a system matrix, $\mathbf{u}$ is the measurement vector, and $\mathbf{R(x)}$ is an (optional) regularization term.
 
# MPIReco comes with a few prepared blueprints, however one can easily add and store their own blueprints. 
# These can be new configurations of existing algorithms and parameters or even new ones developed in other packages. Such packages could provide algorithms which are not based on system-matrix reconstructions. 
# To read more on how to write and integrate such packages, consult the How-Tos.

# ## High-Level Interface
# In the [Getting Started](@ref) section, we created our first MPI reconstruction as follows:
using MPIReco #hide
bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf")) #hide
b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf")) #hide

c = reconstruct("SinglePatch", b;
                   SNRThresh=5,
                   sf = bSF,
                   frames=1:acqNumFrames(b),
                   minFreq=80e3,
                   recChannels=1:rxNumChannels(b),
                   iterations=1,
                   spectralLeakageCorrection=true);

# "SinglePatch" here refers to a blueprint for system-matrix based single-patch image reconstruction. The keyword arguments are used to set parameters defined in the blueprint.

# It is also possible to provide these parameters as a dictionary. 
# This is especially helpful if you want to reuse parameters or change parameters programmatically.
parameters = Dict{Symbol, Any}()
parameters[:SNRThresh] = 5
parameters[:sf] = bSF
parameters[:frames] = 1:acqNumFrames(b)
parameters[:minFreq] = 80e3
parameters[:recChannels] = 1:rxNumChannels(b)
parameters[:iterations] = 1
parameters[:spectralLeakageCorrection] = true;

# You can just pass the dictionary as an argument to the `reconstruct` function as follows:
c2 = reconstruct("SinglePatch", b; parameters...)
isapprox(c.data, c2.data)

# It is also possible to combine both methods of defining parameters:
c3 = reconstruct("SinglePatch", b; parameters..., recChannels = 1:1)
c3.rxNumChannels
# Here, we overwrote the recChannels parameter within in the dictionary.

# We can retrieve a list of all available blueprints that MPIReco can find using:
getRecoPlanList()
# This list returns the filenames of the blueprints. In the reconstruct call, you only need to refer to the name without an extension.


# To get the full paths of each blueprint, use:
planpaths = getRecoPlanList(;full=true)
# You can add custom directories for MPIReco to search for blueprints.
# Other packages might do this automatically when they are loaded.

# It is also possible to specify a direct path to a blueprint in the `reconstruct` function.
index = findfirst(endswith("SinglePatch.toml"), planpaths)
c4 = reconstruct(planpaths[index], b; parameters...)
isapprox(c.data, c4.data)


# ## RecoPlan Interface
# You can also directly access a blueprint via:
plan = MPIRecoPlan("SinglePatch")

# This allows you to see all parameters associated with the blueprint 
# and configure them individually. The keyword arguments previously mentioned 
# apply to all parameters with the same name. Note that different parameters 
# in different contexts can share the same name in RecoPlans.

# To modify specific parameters of the plan, you can directly access its properties:
plan.parameter.reco.sf = bSF

# Although we only updated the system matrix parameter, the plan now also includes 
# derived values from the system matrix, such as the grid size.
plan.parameter.reco

# This worked because the `SinglePatch` blueprint has connections 
# between parameters. These connections depend on the blueprint itself and are not hardcoded in MPIReco.

# You can also utilize the previously created dictionary to configure the plan:
setAll!(plan, parameters)

# Alternatively, you can set all parameters of a specific name directly:
setAll!(plan, :iterations, 3)

# Once the plan is configured, you can construct an algorithm from it using:
algo = build(plan);

# After constructing the algorithm, you can perform multiple reconstructions and reuse the algorithm with:
c = reconstruct(algo, b);
# Depending on the algorithm, it might be faster to reconstruct multiple measurements from a `RecoPlan`, instead of the high-level interface.
# For more information on that, we refer to the How-To for caching.
# Algorithms are usually thread-safe, though they might not necessarily run concurrently.