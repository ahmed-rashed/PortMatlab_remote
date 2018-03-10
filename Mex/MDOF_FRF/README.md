This directory includes Visual Studio Fortran and C++ projects that create Matlab Mex files (shared libraries, i.e.: dll's for Windows). These projects show how you can write code similar to Matlab code but using C++ and Fortran.
Projects in this repositories use:
1) header files within the Mex directory
2) Intel MKL

Till now, projects in this directory are tested on MS Visual Studio 2017, Matlab 2017b and Intel Parallel Studio XE 2017.

Projects under this directory calculate the FRF of a MDOF system using two methods. Details behind these methods are explained in the graduate course AER 631: Dynamics of Structures.
Slow method:
	This includes the following projects:
	1) MDOF_FRF_Visc_Slow_cpp
	2) MDOF_FRF_Visc_Slow_Fortran
	
Fast method:
	This includes the following projects:
	1) MDOF_Eig_Visc_cpp_Armadillo
	2) MDOF_FRF_Visc_cpp_Armadillo
	3) MDOF_Eig_Visc_cpp
	4) MDOF_FRF_Visc_cpp
	5) MDOF_Eig_Visc_Fortran
	6) MDOF_FRF_Visc_Fortran