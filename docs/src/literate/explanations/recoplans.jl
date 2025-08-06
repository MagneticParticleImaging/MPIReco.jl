# # MPIRecoPlans
using MPIReco #hide
# In this explanation we focus on the blueprint mechanism `RecoPlans` of AbstractImageReconstruction and how these are used in MPIReco.
# We start with a quick recap of the RecoPlan interface and then see how MPIReco uses it.

# ## RecoPlans
# `RecoPlans` are blueprints from which parameters and algorithms can be constructed. They are essentially thin-wrappers around nested-key value pairs.
# Plans are able to represent the nested tree structure of algorithms and their parameters. A plan can either be empty, partially or fully parameterized.
# Empty plans are helpful to describe just the structure of an algorithm, which a user can then parameterize.
# If a plan is constructed this way, it is missing all parameters:
using MPIReco.AbstractImageReconstruction
plan = RecoPlan(CommonPreProcessingParameters)
plan.frames

# We can interact with the plan as if it were a mutable variant of the given parameters:
plan.frames = 1:10
plan.bgParams = NoBackgroundCorrectionParameters()
plan

# And we can construct instances of parameters and algorithms from a plan:
parameter = build(plan)

# Likewise we can go from an instance to a plan via:
plan = toPlan(parameter)

# We can either manually empty a plan from its parameters:
plan.frames = missing

# Or use the `clear!` method, which preserves the structure of a blueprint by default:
clear!(plan)

# `RecoPlans` can be written to and read from files:
plan.frames = 1:42
plan.loadasreal = true
plan.numAverages = rand(1:100)
toTOML(stdout, plan) #hide
toTOML(joinpath(@__DIR__, "Parameters.toml"), plan)

# This way one can either serialize a completely parametrized image reconstruction or just the blueprint of an algorithm.

# For more information on RecoPlans, we refer to the How-to section of the AbstractImageReconstruction documentation.

# ## MPIRecoPlan
# MPIReco facilitates access to `RecoPlans` via two additional features.
# The first features deals with loading `RecoPlans`. To load a `RecoPlan` from a file one has to provide the filename and a list of modules in which the respective algorithms and parameters are defined.
# MPIReco tracks a list of directories and modules to grant a user easier access to plans contained within the directories.

# As was shown in the [Implement Reconstruction Packages](@ref) how-to, it is possible to register new modules and directory to MPIReco:
addRecoPlanPath(@__DIR__())
plan = MPIRecoPlan(joinpath(@__DIR__, "Parameters"))

# The second feature was shown in the [Enable Caching](@ref) how-to. This feature is realized by MPIReco keeping a least-recently-used (LRU) cache for recently opened plans.
# As long as the plan structure was not changed in the file or the parameters, a plan can be reused from the cache. While this also saves costs in loading a plan, the main benefit is reusing the same plan instance.
# If the plan in question contains processing steps that implement caching, all algorithms derived from a plan can access this cache. This enables an algorithm to reuse the results of previous processing steps.