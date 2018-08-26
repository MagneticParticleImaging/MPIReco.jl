# Installation and Usage

## Installation

The MPIReco project is a package written in the programming language Julia.
In order to install MPIReco you first have to install [Julia](http://julialang.org/downloads/) in version 0.7 (version 1.0 not yet supported). Then open Julia and enter `]` to open
the package mode. Then enter
```julia
dev https://github.com/MagneticParticleImaging/MPIFiles.jl.git
dev https://github.com/tknopp/LinearSolver.jl.git
dev https://github.com/MagneticParticleImaging/MPIReco.jl.git
```
which will install the package including its dependencies [MPIFiles](https://github.com/MagneticParticleImaging/MPIFiles.jl.git) and [LinearSolver](https://github.com/tknopp/LinearSolver.jl.git). Then enter
```julia
using MPIReco
```
to load the package

In order to obtain the example datasets you have to execute the unit tests which can be done by entering
```julia
test MPIReco
```
within the package mode of Julia.

## Getting Started

In order to get started you can execute
```julia
using Pkg
include(joinpath(dirname(pathof(MPIReco)),"..","test","MultiPatch.jl"))
```
which will perform exemplary reconstructions and plot the results.
