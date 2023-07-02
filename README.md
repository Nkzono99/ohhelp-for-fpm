# ohhelp-for-fpm
OhHelp Library Package: Version 1.1.1 + α for Fortran Package Manager.

## About OhHelp Library Package

OhHelp library is a dynamic load balancing and scalability library that supports massively parallel particle-in-cell simulations using MPI.
Each process handles the particle calculations within the partitioned area while taking on some of the particles handled by other processes so that the overall load is balanced.
Developed by Nakashima et al. [^nakashima_2019] and used mainly in the field of plasma particle simulation, it has been shown to be effective for models such as plasma-satellite interaction and magnetospheric plasma, where the number of particles is time-varying and non-uniform.

The site that provided OhHelp is no longer available, so we are redistributing it and fixing bugs under the license provided by original OhHelp library.

## Install
``` toml
[dependencies]
ohhelp = { git = "https://github.com/Nkzono99/ohhelp-for-fpm" }
```

## References
[^nakashima_2019]: Nakashima, Hiroshi, Yohei Miyake, Hideyuki Usui, and Yoshiharu Omura. 2009. “OhHelp: A Scalable Domain-Decomposing Dynamic Load Balancing for Particle-in-Cell Simulations.” In Proceedings of the 23rd International Conference on Supercomputing, 2009, Yorktown Heights, NY, USA, June 8-12, 2009, 90–99. unknown.
