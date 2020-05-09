Contents of this directory
==========================
  1) header files that enables writing C++ and Fortran codes as easy as Matlab code
	Additionally, these header files eases creation of Matlab Mex files (shared libraries, i.e.: dll's for Windows).
	These header files were tested using MS Visual Studio 2019, Matlab 2019b and Intel Parallel Studio XE 2017.
  2) SDOF_FreeResponse visual studio solution
  3) MDOF_FRF visual studio solution

Instructions to compile and link MEX-Files using Microsoft Visual Studio
========================================================================
1) Create a module definition file as follows
	Project > Add New Item > Module-Definition File (.def)
		Specify the name of the file as "ModuleDefinitionFile.def" (or any other name you like)
		Write the following in the file:
			EXPORTS mexFunction			# This is for a C MEX-file
			or
			EXPORTS _MEXFUNCTION		# This is for a Fortran MEX-file

2) Add your (.c or .cpp or .f90) source files as follows:
	Project > Add Existing Item

3) Modify the project properties as follows:
	Project > Properties > All Configurations >
			General >
				Target Extension > .mexw32 (for win32) or .mexw64 (for win64)
				Configuration Type > Dynamic Library (.dll)

			C/C++ > or Fortran >
				General > Additional Include Directories >
					C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2017.5.267\windows\mkl\include
					C:\Program Files\MATLAB\R2017b\extern\include
					Any other directories you use, e.g.: my Mex folder
					If you are using Fortran 95 MKL functions, add C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2017.5.267\windows\mkl\include\intel64\lp64

				Preprocessor > Preprocessor Definitions >  add MATLAB_MEX_FILE

			Linker >
				General >
					Additional Library Directories: matlabroot\extern\lib\win32\microsoft (for win32) or matlabroot\extern\lib\win64\microsoft (for win64)

				Input >
					Additional Dependencies > libmx.lib, libmex.lib, and libmat.lib, mkl_rt.lib
					Module Definition File > ModuleDefinitionFile.def	(the filename module definition (.def) file you created)

				Debug using C:\Program Files\MATLAB\R2017b\bin\matlab.exe

4) If you are using Intel MKL, follow the instructions for linking with MKL
	These instructions do change from version to another.
	As of March 2018, instructions are available at:
		https://software.intel.com/en-us/articles/intel-math-kernel-library-intel-mkl-compiling-and-linking-with-microsoft-visual-cc
		and
		https://software.intel.com/en-us/mkl-windows-developer-guide-automatically-linking-your-intel-visual-fortran-project-with-intel-mkl
