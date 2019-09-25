# Motion Compensation Example

This folder contains example code for applying motion compensation to multi-patch
MPI data.

## Installation

In order to use this code one first has to download Julia (version 1.x) and install
`MPIReco` by executing

```julia
using Pkg
Pkg.add("MPIReco")
```

You can the switch to the directory of this example by entering
```julia
dir = joinpath(dirname(pathof(MPIReco)), "..","examples","MotionCompensation")
cd(dir)
```

## Execution
After installation and switching to the example directory, the example code can be
executed by entering

```julia
include("example.jl")
```

This will first download all data (about 2 GB) and then perform a reconstruction.
Parameters of the reconstruction are documented in the Julia script and can be
changed. After the reconstruction is done, the script will open a PyPlot window
and show the reconstruction result. Additionally, an animated gif will be stored
into the `img/` folder.
