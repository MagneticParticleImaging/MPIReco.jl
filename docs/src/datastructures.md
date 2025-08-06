# Data Structures

MPIReco contains several different groups of data structures which will be explained in more details in this and the following explanation pages. This page focuses on the core data structures used during image reconstruction, while [MPIRecoPlan](@ref) focuses on the blueprint mechanisms of AbstractImageReconstruction and its interaction with MPIReco. And lastly, [Imaging Operators](@ref) highlights the requirements and structure of the operators used during image reconstruction.

## AbstractImageReconstruction
Image reconstruction using AbstractImageReconstruction.jl allows for flexible control of the reconstruction process and data flow based on individual processing steps. The flow of these steps is defined by Julia's multiple dispatch mechanism applied to algorithms, their (configuration) parameters, and the input arguments of the processing steps.

```julia
abstract type AbstractImageReconstructionParameters end
```
Parameters are adjustable settings used during image reconstruction. Users can typically configure these through keyword arguments. Parameters can be nested and may contain other parameters.

```julia
abstract type AbstractImageReconstructionAlgorithm end
```
Algorithms are the data structures responsible for executing image reconstruction and provide the context within which a given set of parameters is evaluated. Different algorithms may have slightly different implementations of a processing step for the same parameter. Algorithms are usually stateful and implement a thread-safe FIFO behavior.

```julia
function process(algo, params, inputs...)
  # ...
end
```
Processing steps are the internal steps used during image reconstruction. Each algorithm defines the control flow for its processing steps by invoking `process` on itself, a parameter, and some input values, such as a file or an array of measurement data.

A `process` can in turn invoke another `process` function with, for example, a different nested parameter or a changed input value. Processing can occur with either an instance of an algorithm or the type of an algorithm. The latter case represents pure processing steps, which can be easily cached and reused between subsequent image reconstructions.

The [AbstractImageReconstruction.jl documentation](https://juliaimagerecon.github.io/AbstractImageReconstruction.jl) features a complete example of how to assemble all these components into a working image reconstruction algorithm, which provides the high-level interface reconstruct(algo, data). Alternatively, one can read the implementation of the provided MPI reconstruction algorithms.

## MPIReco
MPIReco extends the abstract types of AbstractImageReconstruction with its own type hierarchies for `AbstractMPIRecoAlgorithms <: AbstractReconstructionAlgorithm` and `AbstractMPIRecoParameters <: AbstractReconstructionParameters`. For example, the `reconstruct("SinglePatch", file)` function constructs a `SinglePatchAlgorithm`, which extends `AbstractSinglePatchAlgorithm`, which in turn extends `AbstractMPIRecoAlgorithm`.

The provided parameters are roughly grouped into preprocessing parameters, which focus on loading the correct data from MDFs, and reconstruction parameters, which focus on modifying and using the output of preprocessing to construct and solve the inverse problem of image reconstruction. Finally, post-processing parameters can be applied to the resulting images. While this rough grouping is used in the provided algorithms, custom algorithms can deviate to define and combine parameters as best fits their needs.

Many parameters provide a default implementation `process(algoT::AbstractMPIRecoAlgorithm, param, inputs...)`, allowing any algorithm to invoke a standard utility, such as loading measurement data with a given filtering parameter.
