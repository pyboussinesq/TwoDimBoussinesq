# TwoDimBoussinesq
A pythonic solver for the 2D (x-z) bossinesq equations. It uses streamfunction-vorticity-buoyancy formulation, and solves the equations using pseudo-spectral methods. Time-stepping is performed using a RK3-Theta method.

# pyboussinesq
This 2D code is a first step towards the main goal of the `pyboussinesq` project that aims to build
a pythonic mpi-parallelizable solver for the 3D Boussinesq equations.

# Team
* Gregory L. Wagner
* Cesar B. Rocha

# History
TwoDimBoussinesq started as a script for solving the 2D (x-z) when Boussinesq equations
when Cesar was a fellow of the [WHOI GFD Program](https://www.whoi.edu/gfd/). This script
was crafted in an afternoon hack. Its design partly builds on [pyqg](http://pyqg.readthedocs.org/en/stable/).

The code uses Fourier pseudo-spectral methods. Time-stepping is performed using a semi-implicit, low-storage RK3-theta method.

# Future
The goal of this project is to take this script to the next level, ultimately generalizing the code
to three-dimensions and different basis-functions. The main motivation is to implement the code to study the physics of breaking 
internal gravity waves.

# contributing 
Contribution are very welcome.
