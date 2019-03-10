# Basic Reconstruction

MPIReco.jl provides different reconstruction levels. All of these reconstruction
routines are called `reconstruction` and the dispatch is done based on the input
types.

## On Disk Reconstruction

This is the highest level reconstruction. The function signature is given by
```julia
function reconstruction(d::MDFDatasetStore, study::Study,
                        exp::Experiment, recoParams::Dict)
```
This reconstruction is also called an *on disk* reconstruction because it assumes
that one has a data store (i.e. a structured folder of files) where
the file location is uniquely determined by the study name and experiment number.
All reconstruction parameters are passed to this method by the `recoParams` dictionary.
On disk reconstruction has the advantage that the routine will perform reconstruction
only once for a particular set of parameters. If that parameter set has already
been reconstructed, the data will loaded from disk.
However, the on disk reconstruction needs some experience with dataset stores to
set it up correctly and is not suited for unstructured data.

## In Memory Reconstruction

The next level is the in memory reconstruction. Its function signature reads
```julia
function reconstruction(recoParams::Dict)
```
This routine requires that all parameters are put into a dictionary. An overview
how this dictionary looks like is given in the section [Parameters](@ref).

The above reconstruction method basically does two things
* Pull out the location of measurement data and system matrix from the `recoParams`
  dictionary.
* Pass all parameter to the low level reconstruction method in the form of keyword
  parameters.

In turn the next level reconstruction looks like this
```julia
function reconstruction(bSF::Union{T,Vector{T}}, bMeas::MPIFile; kargs...)
```
There are, however also some reconstruction methods in-between that look like this
```julia
function reconstruction(filenameSF::AbstractString, filenameMeas::AbstractString; kargs...)
function reconstruction(filenameMeas::AbstractString; kargs...)
```
In both cases, an MPIFile is created based on the input filename. The second version
also guesses the system matrix based on what is stored within the measurement
file. This usually only works, if this is executed on a system where the files
are stored at exactly the same location as how they have been measured.

## Middle Layer Reconstruction
