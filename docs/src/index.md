# MPIReco.jl

*Julia package for the reconstruction of magnetic particle imaging (MPI) data*

## Introduction

This project provides algorithms and utility functions for the reconstruction of MPI data. The project
is implemented in the programming language Julia and its algorithms use the interface provided by the [AbstractImageReconstruction.jl](https://github.com/JuliaImageRecon/AbstractImageReconstruction.jl) package.

`MPIReco` contains algorithms for a variety of system matrix based reconstructions:

* Single-Patch Reconstruction
* [Multi-Patch Reconstruction](@ref) for data that has been acquired
  using a focus field sequence
* [Multi-Contrast Reconstruction](@ref) for reconstructions using multiple system matrices
* [Matrix-Compression Techniques](@ref)

Furthermore, the existing algorithms can easily be adapted with new data processing steps and newly created algorithms can be seamlessly integrated into MPIRecos reconstruction interface.

Key features are

* Frequency filtering for memory efficient reconstruction. Only frequencies used
  during reconstructions are loaded into memory.
* Different solvers provided by the package [RegularizedLeastSquares.jl](https://github.com/tknopp/RegularizedLeastSquares.jl)
* High-level and low-level reconstruction interfaces provide maximum flexibility for
  the user
* Flexible algorithm definition and parametrization with [AbstractImageReconstruction.jl](https://github.com/JuliaImageRecon/AbstractImageReconstruction.jl) and the possibilty to define and include custom reconstruction algorithms
* Spectral leakage correction (implemented in
  [MPIFiles.jl](https://github.com/MagneticParticleImaging/MPIFiles.jl))
* Vendor agnostic GPU acceleration

## Installation

Start julia and open the package mode by entering `]`. Then enter
```julia
add MPIReco
```
This will install the packages `MPIReco.jl` and all its dependencies. In particular
this will install the core dependencies [MPIFiles](https://github.com/MagneticParticleImaging/MPIFiles.jl.git), [RegularizedLeastSquares](https://github.com/tknopp/RegularizedLeastSquares.jl.git) and [AbstractImageReconstruction.jl](https://github.com/JuliaImageRecon/AbstractImageReconstruction.jl).

## License / Terms of Usage

The source code of this project is licensed under the MIT license. This implies that
you are free to use, share, and adapt it. However, please give appropriate credit
by citing the project.

## Contact

If you have problems using the software, find mistakes, or have general questions please use
the [issue tracker](https://github.com/MagneticParticleImaging/MPIReco.jl/issues) to contact us.

## Contributors

* [Tobias Knopp](https://www.tuhh.de/ibi/people/tobias-knopp-head-of-institute.html)
* [Martin Möddel](https://www.tuhh.de/ibi/people/martin-moeddel.html)
* [Niklas Hackelberg](https://www.tuhh.de/ibi/people/niklas-hackelberg)
* [Patryk Szwargulski](https://www.tuhh.de/ibi/people/patryk-szwargulski.html)
