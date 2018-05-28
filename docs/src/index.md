# MPI Reco

Julia package for reconstruction of magnetic particle imaging (MPI) data

## Introduction

This project provides functions for the reconstruction of MPI data. The project
is implemented in the programming language Julia and contains algorithms for

* regular 1D / 2D / 3D image reconstruction using a system matrix based approach
* multi-patch and multi-gradient reconstruction for data the has been acquired
  using a focus field sequence
* multi-colored image reconstruction
* matrix-compression techniques

Key features are

* frequency filtering for memory efficient reconstruction. Only frequencies used
  during reconstructions are loaded into memory.
* spectral leakage correction (implemented in MPIFiles.jl)


## License / Terms of Usage

The source code of this project is licensed under the MIT license. This implies that
you are free to use, share, and adapt it. However, you must give appropriate credit by citing the project.


## Contact

If you have problems using the software, find mistakes, or have general questions please use
the [issue tracker](https://github.com/MagneticParticleImaging/MPIReco.jl/issues) to contact us.

## Contributors

* [Florian Griese](https://www.tuhh.de/ibi/people/florian-griese.html)
* [Nadine Gdaniec](https://www.tuhh.de/ibi/people/nadine-gdaniec.html)
* [Tobias Knopp](https://www.tuhh.de/ibi/people/tobias-knopp-head-of-institute.html)
* [Martin Möddel](https://www.tuhh.de/ibi/people/martin-moeddel.html)
* [Patryk Szwargulski](https://www.tuhh.de/ibi/people/patryk-szwargulski.html)
