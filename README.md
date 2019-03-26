# frehdc

FrehdC (Fine Resolution Environmental Hydrodynamic Model in C) is a 2D depth-integrated semi-implicit finite-volume numerical model that simulates free surface flow and scalar transport in shallow estuaries, river deltas and coastal marshes. The original Frehd model was written by Dr. Ben Hodges from the University of Texas at Austin. This 2D C-version Frehd is created by Zhi Li to improve its computation efficiency (by enabling parallel computing with MPI). The solution algorithm follows the study by Casulli (see references below). The linear system is solved with the laspack package by Tomas Skalicky (http://www.mgnet.org/mgnet/Codes/laspack/html/laspack.html).

The current available version of FrehdC (v4.2) is still at testing stage. The comments and documentations may not be user-friendly. For any questions regarding the model please ask Zhi Li at zhili@utexas.edu.


# references:

Casulli V (1990), Semi-implicit finite-difference methods for the 2-dimensional shallow-water equations, Journal of Computational Physics, 86: 56–74.

Casulli V, Cattani E (1994), Stability, accuracy and efficiency of a semi-implicit method for three-dimensional shallow water flow, Computer & Mathematics with Applications, 27(4): 99–112.

Casulli V (1999), A semi-implicit finite difference method for non-hydrostatic free-surface flows,International Journal for Numerical Methods in Fluids, 30: 425–440.

Hodges BR, Imberger J, Saggio A, Winters K (2000), Modeling basin-scale internal waves in a stratified lake, Limnology and Oceanography, 47(7): 1603–1620.

Hodges BR (2004), Accuracy order of Crank-Nicolson discretization for hydrostatic free surface flow, ASCE Journal of Engineering Mechanics, 130(8): 904–910, DOI: 10.1061/(ASCE)0733- 9399(2004)130:8(904)

Hodges BR, Rueda FJ (2008), Semi-implicit two-level predictor-corrector methods for non-linearly coupled, hydrostatic, barotropic/baroclinic flows, International Journal of Computational Fluid Dynamics, 22(9): 593–607, http://dx.doi.org/10.1080/10618560802353389

Hodges BR (2015), Representing hydrodynamically important blocking features in coastal or riverine lidar topography, Natural Hazards & Earth System Sciences 56: 1011–1023, http://dx.doi.org/10.5194/nhessd-3-1427-2015

Rueda FJ, Sanmiguel-Rojas E, Hodges BR (2007),Baroclinic stability for a family of two-level,semi-implicit numerical methods for the 3D shallow water equations, International Journal of Numerical Methods in Fluids, 54(3): 237–268. http://dx.doi.org/10.1002/fld.1391

Wadzuk BM, Hodges BR (2009), Hydrostatic versus nonhydrostatic Euler-Equation modeling of nonlinear internal waves, ASCE Journal of Engineering Mechanics, 135(4), 1069–1080, http://dx.doi.org/10.1061/(ASCE)0733-9399(2009)135:10(1069)
