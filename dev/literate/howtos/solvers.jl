# # Change and Configure Solvers
include("../../download.jl") #hide
using MPIReco #hide
bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf")) #hide
b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf")); #hide
# MPIReco uses [RegularizedLeastSquares.jl](https://github.com/JuliaImageRecon/RegularizedLeastSquares.jl) as its optimization backend.
# As such it can use any of the solvers and regularization terms defined for and in RegularizedLeastSquares.

# In general, RegularizedLeastSquares aims to solve problems of the form:
# ```math
# \begin{equation}
#   \underset{\mathbf{c}}{argmin} \frac{1}{2}\vert\vert \mathbf{S}\mathbf{c}-\mathbf{u} \vert\vert_2^2 + + \mathbf{R(x)}
# \end{equation}
# ```
# where $\mathbf{S}$ is a system matrix operator and $\mathbf{u}$ is the measurement vector. Both are provided by MPIReco and constructed based on the
# specific parameters and processing steps of a given reconstruction algorithm. The (optional) regularization term $\mathbf{R(x)}$ can be used to encode additional prior information about the inverse problem.
# Different solvers of RegularizedLeastSquares can work with different regularization terms. For example, the Kaczmarz and CGNR solver are defined for the $l^2_2$-norm.
# More information on available solvers and regularization terms can be found RegularizedLeastSquares documentation.

# ## Solver and Solver Parameters
# MPIReco offers different ways to configure the used solvers, similar to how it for example offers different weighting strategies.
# If we take a look at the default solver parameters from the single-patch reconstruction:
plan = MPIRecoPlan("SinglePatch")
solverParams = plan.parameter.reco.solverParams
# We see that the default parameters are of type `ElaborateSolverParameters`. These are a special set of parameters, because they contain the union of all solver
# parameters defined in RegularizedLeastSquares. The current selection of parameters changes based on the selected solver:
solverParams.solver = Kaczmarz
solverParams
# or:
solverParams.solver = FISTA
solverParams
# With these parameters we can define both the solver and the solver parameters as defined in RegularizedLeastSquares.

# It is also possible to add new solvers, such as variants of existing solver specialised for MPI.
# To do this, we have to implement a new solver variant for RegularizedLeastSquares. We can either implement a new solver or implement a new variant of a solver via its state:
using MPIReco.RegularizedLeastSquares
mutable struct OurSolver{OP, T} <: RegularizedLeastSquares.AbstractLinearSolver
  S::OP
  c::Vector{T}
  u::Vector{T}
  notification::String
end
function OurSolver(S; notification = "Test", kwargs...)
  u = zeros(eltype(S), size(S, 1))
  c = zeros(eltype(S), size(S, 2))
  return OurSolver(S, c, u, notification)
end
function RegularizedLeastSquares.init!(solver::OurSolver{OP, T}, u) where {OP, T}
  solver.u .= u
  solver.c .= zero(T)
end
function Base.iterate(solver::OurSolver)
  @info solver.notification
  solver.c .= solver.S \ solver.u
  return nothing
end
RegularizedLeastSquares.solversolution(solver::OurSolver) = solver.c

# Afterwards, we can just use the solver type as usual:
params = Dict{Symbol, Any}()
params[:SNRThresh] = 5
params[:frames] = 1:1
params[:minFreq] = 80e3
params[:recChannels] = 1:2
params[:spectralLeakageCorrection] = true
params[:sf] = bSF;

c = reconstruct("SinglePatch", b; params..., solver = OurSolver)
using CairoMakie #hide
fig = heatmap(c[1, :, :, 1, 1].data.data)
hidedecorations!(fig.axis)
fig

# If our custom solver requires custom parameters, we could implement custom solver parameters in MPIReco:
Base.@kwdef struct OurSolverParameters <: MPIReco.AbstractSolverParameters{OurSolver}
  notification::String
  enforceReal::Bool=true
  enforcePositive::Bool=true
end
reconstruct("SinglePatch", b; params..., solverParams = OurSolverParameters(notification = "Custom"));
# More on custom parameters for MPIReco can be found in [Custom Data Processing and Algorithms](@ref).

# ## Regularization Term
# RegularizedLeastSquares allows for the generartion of flexible regularization using a vector of (nested) regularization terms. 
# By providing such a vector, we can configure the regularization term similar to the solver:
setAll!(plan, :reg, [L2Regularization(0.1)])
# Because MPI oftens constrain its solutions to positive, real numbers, MPIReco offers some shortcurts to enable/disable these constraints:
setAll!(plan, :enforceReal, true)
# Internally, these are just additional regularization terms added to the provided regularization terms. For more complex regularization
# examples, take a look at the RegularizedLeastSquares documentation.

# To implement a custom regularization we have to implement a proximal mapping:
# ```math
# \begin{equation}
#   prox_g (\mathbf{x}) = \underset{\mathbf{u}}{argmin} \frac{1}{2}\vert\vert \mathbf{u}-\mathbf{x} \vert {\vert}^2 + g(\mathbf{x}).
# \end{equation}
# ```
# For many regularizers, the proximal map can be computed efficiently in a closed form.

# In order to implement these proximal mappings, we have to defines the following type functions for our new proximal map:
# ```julia
# abstract type AbstractRegularization
# prox!(reg::AbstractRegularization, x)
# norm(reg::AbstractRegularization, x)
# ```
# Here `prox!(reg, x)` is an in-place function which computes the proximal map on the input-vector `x`. The function `norm` computes the value of the corresponding term in the inverse problem. RegularizedLeastSquares.jl provides `AbstractParameterizedRegularization` and `AbstractProjectionRegularization` as core regularization types.
