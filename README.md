# RectCylinderSIMPLERTVD
Solution of the incompressible navier-stokes equations using FVM with a cartesian uniform grid.

It solves a 3D incompressible flow with constant viscosity around a rectangle cylinder.

The main source for this work is from "An Introduction to Computational Fluid Dynamics: The Finite Volume Method" - Veersteeg and Malalasekera.

The pressure-velocity coupling is made using the paper "IMPROVED NUMERICAL CODES FOR SOLVING THREE-DIMENSIONAL UNSTEADY FLOWS" - Chin-Yuan Perng and Robert L. Street

IMPORTANT: TO SOLVE THE SPARSE LINEAR SYSTEMS, I AM USING THE AMGCL (google it), IT'S NOT INCLUDED IN THIS PROJECT. IF YOU WANT TO COMPILE IT, YOU WILL NEED IT. 
