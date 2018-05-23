# Installation and Usage

## Installation

The MPIReco project is a package written in the programming language Julia.
In order to install MPIReco you first have to install [Julia](http://julialang.org/downloads/) in version 0.6. Then open Julia and enter
```julia
Pkg.clone("https://github.com/MagneticParticleImaging/MPIReco.jl.git")
```
which will install the package. Then enter
```julia
using MPIReco
```
which will install the dependencies [MPIFiles](https://github.com/MagneticParticleImaging/MPIFiles.jl.git) and [LinearSolver](https://github.com/tknopp/LinearSolver.jl.git). Further dependencies have already been installed during the clone of the package.

In order to obtain the example datasets you have to execute the unit tests which can be done by entering
```julia
Pkg.test("MPIReco")
```

## Getting Started

In order to get started you can execute on
```julia
include(Pkg.dir("MPIReco","test/MultiPatch.jl"))
```
which will perform exemplary reconstructions and plot the results.
