# This make file assumes the following environment variables are defined
# MATLABROOT

ifeq ($(OS),Windows_NT)
	Ext=mexw64
	ObjExt=obj
	PWD :=$(notdir $(shell cd))
	mkdir_CMD=mkdir
	rm_CMD=rmdir /S /Q
else
	Ext=mexa64
	ObjExt=o
	PWD :=$(notdir $(shell pwd))
	mkdir_CMD=mkdir -p
	rm_CMD=rm -rf
endif

EXEC_main =$(PWD).$(Ext)

SRC_C=$(wildcard *.c)
SRC_CPP=$(wildcard *.cpp)
SRC_F=$(wildcard *.f90)

MACROS=-DMX_COMPAT_64 -DMATLAB_MEXCMD_RELEASE=R2018a -DUSE_MEX_CMD -DMATLAB_MEX_FILE

#Note that $(OutDebugDir): and $(OutReleaseDir): are targets. Therefore, make sure that any of them does not conflict with other targets, especially debug: and release:
OutRootDir=.
DebugDir =Debug
ReleaseDir =Release

OBJ_C_Debug   =$(patsubst %.c,  $(OutRootDir)/$(DebugDir)/%.$(ObjExt),$(SRC_C))
OBJ_CPP_Debug =$(patsubst %.cpp,$(OutRootDir)/$(DebugDir)/%.$(ObjExt),$(SRC_CPP))
OBJ_F_Debug   =$(patsubst %.f90,$(OutRootDir)/$(DebugDir)/%.$(ObjExt),$(SRC_F))
Mod_F_Debug   =$(OutRootDir)/$(DebugDir)/mx_interfaces.mod $(OutRootDir)/$(DebugDir)/mex_interfaces.mod

OBJ_C_Release   =$(patsubst %.c,  $(OutRootDir)/$(ReleaseDir)/%.$(ObjExt),$(SRC_C))
OBJ_CPP_Release =$(patsubst %.cpp,$(OutRootDir)/$(ReleaseDir)/%.$(ObjExt),$(SRC_CPP))
OBJ_F_Release   =$(patsubst %.f90,$(OutRootDir)/$(ReleaseDir)/%.$(ObjExt),$(SRC_F))
Mod_F_Release   =$(OutRootDir)/$(ReleaseDir)/mx_interfaces.mod $(OutRootDir)/$(ReleaseDir)/mex_interfaces.mod

# Do not delete the $(Mod_F_Release) $(Mod_F_Debug) files
.SECONDARY: $(Mod_F_Release) $(Mod_F_Debug)

.PHONY: MSVC_IVF_Release MSVC_IVF_Debug GCC_Release GCC_Debug MINGW_Release MINGW_Debug cleanDebug cleanRelease clean

CC=gcc
CXX=g++
Compile_Flag=-c
CompileOutFlag=-o 
FC=gfortran
LD=gcc
CFLAGS =-I"$(MATLABROOT)/extern/include"
CPPFLAGS =$(CFLAGS)
FFLAGS =$(CFLAGS)
LDFLAGS  =-Wl,--no-undefined -shared
LinkOutFlag =-o 
ModOutFlag =-J

################################################
MSVC_IVF_CC:=cl
MSVC_IVF_CXX:=cl
MSVC_IVF_FC:=ifort
MSVC_IVF_LD:=link
MSVC_IVF_MACROS:=$(patsubst -D%,/D%,$(MACROS)) /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0
MSVC_IVF_CFLAGS   :=$(CFLAGS) /Zp8 /GR /W3 /EHs /nologo /MD
MSVC_IVF_CPPFLAGS :=$(CPPFLAGS) /Zp8 /GR /W3 /EHs /nologo /MD
MSVC_IVF_FFLAGS   :=$(FFLAGS) /nologo /fpp /Qprec /MD /fp:source /assume:bscc /integer-size:64
MSVC_IVF_C_LIBS :=libmx.lib libmex.lib libmat.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib
MSVC_IVF_CPP_LIBS :=libMatlabDataArray.lib libMatlabEngine.lib
MSVC_IVF_F_LIBS :=
MSVC_IVF_LibPath :=/LIBPATH:"$(MATLABROOT)/extern/lib/win64/microsoft"
MSVC_IVF_LDFLAGS :=/nologo /manifest /DLL /DEF:ModuleDefinitionFile.def
MSVC_IVF_LinkOutFlag :=/out:
MSVC_IVF_ModOutFlag :=/module:
MSVC_IVF_Compile_Flag:=/c
MSVC_IVF_CompileOutFlag:=/Fo

MSVC_IVF_Release: CC:=$(MSVC_IVF_CC)
MSVC_IVF_Release: CXX:=$(MSVC_IVF_CXX)
MSVC_IVF_Release: FC:=$(MSVC_IVF_FC)
MSVC_IVF_Release: LD:=$(MSVC_IVF_LD)
MSVC_IVF_Release: C_LIBS:=$(MSVC_IVF_C_LIBS)
MSVC_IVF_Release: CPP_LIBS:=$(MSVC_IVF_CPP_LIBS)
MSVC_IVF_Release: F_LIBS:=$(MSVC_IVF_F_LIBS)
MSVC_IVF_Release: LibPath:=$(MSVC_IVF_LibPath)
MSVC_IVF_Release: LDFLAGS:=$(MSVC_IVF_LDFLAGS)
MSVC_IVF_Release: LinkOutFlag:=$(MSVC_IVF_LinkOutFlag)
MSVC_IVF_Release: ModOutFlag:=$(MSVC_IVF_ModOutFlag)
MSVC_IVF_Release: MACROS   :=$(MSVC_IVF_MACROS) /DNDEBUG
MSVC_IVF_Release: CFLAGS   :=$(MSVC_IVF_CFLAGS) /O2 /Oy-
MSVC_IVF_Release: CPPFLAGS :=$(MSVC_IVF_CPPFLAGS) /O2 /Oy-
MSVC_IVF_Release: FFLAGS   :=$(MSVC_IVF_FFLAGS) /O2
MSVC_IVF_Release: Compile_Flag:=$(MSVC_IVF_Compile_Flag)
MSVC_IVF_Release: CompileOutFlag:=$(MSVC_IVF_CompileOutFlag)

MSVC_IVF_Debug: CC:=$(MSVC_IVF_CC)
MSVC_IVF_Debug: CXX:=$(MSVC_IVF_CXX)
MSVC_IVF_Debug: FC:=$(MSVC_IVF_FC)
MSVC_IVF_Debug: LD:=$(MSVC_IVF_LD)
MSVC_IVF_Debug: C_LIBS:=$(MSVC_IVF_C_LIBS)
MSVC_IVF_Debug: CPP_LIBS:=$(MSVC_IVF_CPP_LIBS)
MSVC_IVF_Debug: F_LIBS:=$(MSVC_IVF_F_LIBS)
MSVC_IVF_Debug: LibPath:=$(MSVC_IVF_LibPath)
MSVC_IVF_Debug: LDFLAGS  :=$(MSVC_IVF_LDFLAGS) /debug /PDB:"$(OutRootDir)/$(DebugDir)/$(EXEC_main).mexw64.pdb"
MSVC_IVF_Debug: LinkOutFlag:=$(MSVC_IVF_LinkOutFlag)
MSVC_IVF_Debug: ModOutFlag:=$(MSVC_IVF_ModOutFlag)
MSVC_IVF_Debug: MACROS   :=$(MSVC_IVF_MACROS)
MSVC_IVF_Debug: CFLAGS   :=$(MSVC_IVF_CFLAGS) /Z7
MSVC_IVF_Debug: CPPFLAGS :=$(MSVC_IVF_CPPFLAGS) /Z7
MSVC_IVF_Debug: FFLAGS   :=$(MSVC_IVF_FFLAGS) /Z7
MSVC_IVF_Debug: Compile_Flag:=$(MSVC_IVF_Compile_Flag)
MSVC_IVF_Debug: CompileOutFlag:=$(MSVC_IVF_CompileOutFlag)
################################################

################################################
GCC_MACROS   :=$(MACROS) -D_GNU_SOURCE
GCC_CFLAGS   :=$(CFLAGS) -fexceptions -fno-omit-frame-pointer -fPIC -pthread
GCC_CPPFLAGS :=$(GCC_CFLAGS)
GCC_FFLAGS   :=$(FFLAGS) -cpp -fPIC -fbackslash -fdefault-integer-8
GCC_base_LIBS :=-lmx -lmex -lmat -lm
GCC_C_LIBS :=$(GCC_base_LIBS) -lstdc++
GCC_CPP_LIBS :=-lMatlabDataArray
GCC_F_LIBS :=-lgfortran
GCC_LibPath :=-L"${MATLABROOT}/bin/glnxa64"
GCC_LDFLAGS :=$(LDFLAGS) -pthread -Wl,--version-script,"mexFunction.map" -Wl,-rpath-link,"${MATLABROOT}/bin/glnxa64"

GCC_Release: C_LIBS:=$(GCC_C_LIBS)
GCC_Release: CPP_LIBS:=$(GCC_CPP_LIBS)
GCC_Release: F_LIBS:=$(GCC_F_LIBS)
GCC_Release: LibPath:=$(GCC_LibPath)
GCC_Release: MACROS   :=$(GCC_MACROS) -DNDEBUG
GCC_Release: CFLAGS   :=$(GCC_CFLAGS) -O2 -fwrapv
GCC_Release: CPPFLAGS :=$(GCC_CPPFLAGS) -O2 -fwrapv
GCC_Release: FFLAGS   :=$(GCC_FFLAGS) -O2
GCC_Release: LDFLAGS  :=$(GCC_LDFLAGS) -O

GCC_Debug: MACROS:=$(GCC_MACROS)
GCC_Debug: C_LIBS:=$(GCC_C_LIBS)
GCC_Debug: CPP_LIBS:=$(GCC_CPP_LIBS)
GCC_Debug: F_LIBS:=$(GCC_F_LIBS)
GCC_Debug: LibPath:=$(GCC_LibPath)
GCC_Debug: CFLAGS   :=$(GCC_CFLAGS) -g
GCC_Debug: CPPFLAGS :=$(GCC_CPPFLAGS) -g
GCC_Debug: FFLAGS   :=$(GCC_FFLAGS) -g
GCC_Debug: LDFLAGS  :=$(GCC_LDFLAGS) -g
################################################

################################################
MINGW_CFLAGS   :=$(CFLAGS) -fexceptions -fno-omit-frame-pointer -m64
MINGW_CPPFLAGS :=$(MINGW_CFLAGS)
MINGW_C_LIBS :=-llibmx -llibmex -llibmat
MINGW_CPP_LIBS :=-lm -llibMatlabDataArray -llibMatlabEngine
MINGW_LibPath :=-L"${MATLABROOT}/extern/lib/win64/mingw64"
MINGW_LDFLAGS :=$(LDFLAGS) -m64 -Wl,"ModuleDefinitionFile.def"

MINGW_Release: C_LIBS:=$(MINGW_C_LIBS)
MINGW_Release: CPP_LIBS:=$(MINGW_CPP_LIBS)
MINGW_Release: LibPath:=$(MINGW_LibPath)
MINGW_Release: MACROS   :=$(MINGW_MACROS) -DNDEBUG
MINGW_Release: CFLAGS   :=$(MINGW_CFLAGS) -O2 -fwrapv
MINGW_Release: CPPFLAGS :=$(MINGW_CPPFLAGS) -O2 -fwrapv
MINGW_Release: LDFLAGS  :=$(MINGW_LDFLAGS) -s

MINGW_Debug: C_LIBS:=$(MINGW_C_LIBS)
MINGW_Debug: CPP_LIBS:=$(MINGW_CPP_LIBS)
MINGW_Debug: LibPath:=$(MINGW_LibPath)
MINGW_Debug: CFLAGS   :=$(MINGW_CFLAGS) -g
MINGW_Debug: CPPFLAGS :=$(MINGW_CPPFLAGS) -g
MINGW_Debug: LDFLAGS  :=$(MINGW_LDFLAGS) -g
################################################

################################################
#Callable Targets
MSVC_IVF_Release GCC_Release MINGW_Release: $(OutRootDir)/$(ReleaseDir)/$(EXEC_main)

MSVC_IVF_Debug GCC_Debug MINGW_Debug:       $(OutRootDir)/$(DebugDir)/$(EXEC_main)

cleanDebug:
	-$(rm_CMD) "$(OutRootDir)/$(DebugDir)"

cleanRelease:
	-$(rm_CMD) "$(OutRootDir)/$(ReleaseDir)"

clean: cleanDebug cleanRelease
############################################################
# Pattern Rules
$(OutRootDir)/$(ReleaseDir)/%.$(ObjExt): %.c | $(OutRootDir)/$(ReleaseDir)
	$(CC) $(Compile_Flag) $(MACROS) $(CFLAGS) $< $(CompileOutFlag)"$@"

$(OutRootDir)/$(DebugDir)/%.$(ObjExt):   %.c | $(OutRootDir)/$(DebugDir)
	$(CC) $(Compile_Flag) $(MACROS) $(CFLAGS) $< $(CompileOutFlag)"$@"

$(OutRootDir)/$(ReleaseDir)/%.$(ObjExt): %.cpp | $(OutRootDir)/$(ReleaseDir)
	$(CXX) $(Compile_Flag) $(MACROS) $(CPPFLAGS) $< $(CompileOutFlag)"$@"

$(OutRootDir)/$(DebugDir)/%.$(ObjExt): %.cpp | $(OutRootDir)/$(DebugDir)
	$(CXX) $(Compile_Flag) $(MACROS) $(CPPFLAGS) $< $(CompileOutFlag)"$@"

$(OutRootDir)/$(ReleaseDir)/%.$(ObjExt): %.f90 | $(Mod_F_Release) $(OutRootDir)/$(ReleaseDir)
	$(FC) $(Compile_Flag) $(MACROS) $(FFLAGS) $< $(ModOutFlag)"$(OutRootDir)/$(ReleaseDir)" $(CompileOutFlag)"$@"

$(OutRootDir)/$(DebugDir)/%.$(ObjExt): %.f90 | $(Mod_F_Debug) $(OutRootDir)/$(DebugDir)
	$(FC) $(Compile_Flag) $(MACROS) $(FFLAGS) $< $(ModOutFlag)"$(OutRootDir)/$(DebugDir)" $(CompileOutFlag)"$@"

$(OutRootDir)/$(ReleaseDir)/%.mod: ../../%.f90 | $(OutRootDir)/$(ReleaseDir)
	$(FC) $(Compile_Flag) $(MACROS) $(FFLAGS) $< $(ModOutFlag)"$(OutRootDir)/$(ReleaseDir)" $(CompileOutFlag)"$(OutRootDir)/$(ReleaseDir)/$*.$(ObjExt)"

$(OutRootDir)/$(DebugDir)/%.mod: ../../%.f90 | $(OutRootDir)/$(DebugDir)
	$(FC) $(Compile_Flag) $(MACROS) $(FFLAGS) $< $(ModOutFlag)"$(OutRootDir)/$(DebugDir)" $(CompileOutFlag)"$(OutRootDir)/$(DebugDir)/$*.$(ObjExt)"

############################################################

$(OutRootDir)/$(ReleaseDir)/$(EXEC_main): $(OBJ_C_Release) $(OBJ_CPP_Release) $(OBJ_F_Release)
	$(LD) $(LDFLAGS) $(OBJ_C_Release) $(OBJ_CPP_Release) $(OBJ_F_Release) $(LibPath) $(C_LIBS) $(CPP_LIBS) $(F_LIBS) $(LinkOutFlag)$@

$(OutRootDir)/$(DebugDir)/$(EXEC_main): $(OBJ_C_Debug) $(OBJ_CPP_Debug) $(OBJ_F_Debug)
	$(LD) $(LDFLAGS) $(OBJ_C_Debug) $(OBJ_CPP_Debug) $(OBJ_F_Debug) $(LibPath) $(C_LIBS) $(CPP_LIBS) $(F_LIBS) $(LinkOutFlag)$@

$(OutRootDir)/$(ReleaseDir) $(OutRootDir)/$(DebugDir):
	-$(mkdir_CMD) $@
