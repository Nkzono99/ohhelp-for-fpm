# ohhelp-for-fpm
OhHelp Library Package: Version 1.1.1 + α for Fortran Package Manager.

## About OhHelp Library Package

OhHelp library is a dynamic load balancing and scalability library that supports massively parallel particle-in-cell simulations using MPI.
Each process handles the particle calculations within the partitioned area while taking on some of the particles handled by other processes so that the overall load is balanced.
Developed by Nakashima et al. [^1] and used mainly in the field of plasma particle simulation, it has been shown to be effective for models such as plasma-satellite interaction and magnetospheric plasma, where the number of particles is time-varying and non-uniform.

The site that provided OhHelp is no longer available, so we are redistributing it and fixing bugs under the license provided by original OhHelp library.

## Install
``` toml
[dependencies]
ohhelp = { git = "https://github.com/Nkzono99/ohhelp-for-fpm" }
```

## References
[^1] Nakashima, Hiroshi. n.d. “OhHelp Library Package for Scalable Domain-Decomposed PIC Simulation.”
