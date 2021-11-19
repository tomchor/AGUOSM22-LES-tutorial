# AGUOSM22-LES-tutorial

Welcome to the online resource for the LES tutorial at the 2022 AGU Ocean Sciences Meeting. 


## Abstract and intro 

Large eddy simulations (LES) are an increasingly important tool for studying turbulent ocean surface and bottom boundary layers. LES approximate the smallest scales of motion with an analytic turbulence closure, thereby dramatically reducing its computational cost compared to direct numerical simulation (where all scales of motion are resolved) while producing a solution with much higher fidelity than “Reynolds Averaged” turbulence models used in general circulation modeling.

LES has been applied to oceanic problems for decades, but its application has been limited by limited software availability and computational costs. Today, the increased availability of LES software, faster CPUs, increased availability of GPUs and GPU-capable software, and increasing scientific interest in the effects of small-scale turbulence of the ocean and Earth system mean that LES will likely become increasingly common in physical oceanography.


## Outline

In this tutorial we will:

- review basic principles of LES
- basics of ocean LES
  - closures used in ocean simulations
  - boundary conditions
  - sponge layers
  - models for the effects of surface waves
- review common pitfalls and best practices
  - dependence of subgrid-scale fluxes on resolution
  - spurious dynamics



## Resources

Resources for learning the theory of large-eddy simulations:

- [Scholarpedia article about subgrid-scale models](http://www.scholarpedia.org/article/Turbulence:_Subgrid-Scale_Modeling) curated by Prof. Charles Meneveau
- Chapter 13 of [Pope's book](https://www.cambridge.org/highereducation/books/turbulent-flows/C58EFF59AF9B81AE6CFAC9ED16486B3A#overview)
- Deardorffs's [1970 paper on turblent channel flows](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/numerical-study-of-threedimensional-turbulent-channel-flow-at-large-reynolds-numbers/D84769F4A3443E4C87E8878303890999)
- Section 3 of Peter Sullivan and James McWilliam's [paper "Dynamics of Winds and Currents Coupled to Surface Waves"](https://www-annualreviews-org.proxy-um.researchport.umd.edu/doi/abs/10.1146/annurev-fluid-121108-145541).
- Section 2.1 of [Chamecki et al.'s review paper on material transport in the ocean surface boundary layer](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019RG000655).


Some LES-capable ocean-friendly open-source software:

- [Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/) 
  - very flexible user-friendly Julia code that runs on GPUs
- [PALM](https://palm.muk.uni-hannover.de/trac)
  - code for running atmospheric and oceanic boundary-layer flows in massively parallel CPU architectures
- [DIABLO](https://github.com/johnryantaylor/DIABLO) and [LESGO](https://lesgo.me.jhu.edu/)
  - both are pseudo-spectral codes writen in Fortran with LES closures for doubly periodic domains


Some LES-capable open-source software that isn't geared towards ocean simulations. These software may have a steeper learning curve, but are more flexible:

- [OpenFOAM](https://www.openfoam.com/)
  - flexible code written in C++ with a very active user base
- [Fluidity](https://fluidityproject.github.io/)
  - flexible code with unstructured adaptve mesh capabilities and a set of LES closures to choose from
