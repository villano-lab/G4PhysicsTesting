# --------------------------------------------------------------
# GNUmakefile for physics list user.  
# JPW. Fri Jul 25 10:39:58 CEST 2003
# --------------------------------------------------------------
TOPDIR=../
SIMEVENTDIR=$(TOPDIR)getSimEvent/

name := NeutReflectometry

G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

#use the following line to remove compilation of visualization
#G4VIS_NONE := 1 

MYFILELIBS   = -L$(SIMEVENTDIR)  -lgetSimEvent 
MYFILEFLAGS  = -I$(SIMEVENTDIR)

#--- The matlab libraries ....

#MATLABLIBS   = -L/Xapps/Number_Crunching/Matlab_r2007b/extern/lib/maci/  -lmx -lut -lmat
#MORELIBS = -L/Xapps/Number_Crunching/Matlab_r2007b/bin/maci/ -lut -lmat -lmx
#ANOTHERONE = -L/Xapps/Number_Crunching/Matlab_r2007b/sys/os/maci -ldl -lXm
#MATLABFLAGS  = -I/Xapps/Number_Crunching/Matlab_r2007b/extern/include

CPPFLAGS += `root-config --cflags`
EXTRALIBS += `root-config --libs`
#CPPFLAGS += $(MYFILEFLAGS)
#EXTRALIBS += $(MYFILELIBS) 
#CPPFLAGS += $(MATLABFLAGS)
#EXTRALIBS += $(MATLABLIBS) $(ANOTHERONE) $(MORELIBS)
#EXTRALIBS += -licudata.24
#EXTRALIBS += -licui18n.24
#EXTRALIBS += -licuuc.24
#EXTRALIBS += -lustdio.24
#EXTRALIBS += -lz
ifdef G4EXPAT_PATH
  EXTRALIBS += -L$(G4EXPAT_PATH)
endif

include $(G4INSTALL)/config/architecture.gmk

#
# define G4LISTS_BASE, if you have your own physics lists area installed point
# G4LISTS_BASE to the directory, that contains the subdirectory 'hadronic'.
#
#ifndef G4LISTS_BASE
#  EXTRALIBS += -L$(G4LIB)/plists/$(G4SYSTEM)
#  G4LISTS_BASE = $(G4INSTALL)/physics_lists
#else
#  EXTRALIBS += -L$(G4LISTS_BASE)/hadronic/plists/lib/$(G4SYSTEM)
#endif
#EXTRALIBS += -lG4parameterisation
#
# Select your physics lists to link against.
#
#EXTRALIBS += -lLBE
##EXTRALIBS += -lG4shortlived

#CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/LBE/include
#CPPFLAGS += -I$(G4LISTS_BASE)/hadronic/Packaging/include

.PHONY: all
 
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

