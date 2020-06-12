Instructions to compile and link MEX-Files using Microsoft Visual Studio
========================================================================
Reference [www.mathworks.com/help/releases/R2017b/matlab/matlab_external/compiling-mex-files-with-the-microsoft-visual-c-ide.html]
This assumes you have initially defined the following environment variables:
    MATLABROOT

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
					$(MKLIncludeDir)
					$(MATLABROOT)/extern/include

				Preprocessor > Preprocessor Definitions >  add:
                    MX_COMPAT_64;MATLAB_MEXCMD_RELEASE=R2018a;USE_MEX_CMD;MATLAB_MEX_FILE

			Linker >
				General >
					Additional Library Directories: $(MATLABROOT)/extern/lib/win32/microsoft (for win32) or $(MATLABROOT)/extern/lib/win64/microsoft (for win64)

				Input >
					Additional Dependencies > libmx.lib;libmex.lib;libmat.lib;libMatlabDataArray.lib;libMatlabEngine.lib
					Module Definition File > ModuleDefinitionFile.def	(the filename module definition (.def) file you created)

				Debug using $(MATLABROOT)/bin/matlab.exe

4) If you are using Intel MKL, follow the instructions for linking with MKL
	These instructions do change from version to another.
	As of March 2018, instructions are available at:
		https://software.intel.com/en-us/articles/intel-math-kernel-library-intel-mkl-compiling-and-linking-with-microsoft-visual-cc
		and
		https://software.intel.com/en-us/mkl-windows-developer-guide-automatically-linking-your-intel-visual-fortran-project-with-intel-mkl
