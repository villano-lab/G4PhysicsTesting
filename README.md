# G4PhysicsTesting -- Small Geant4 simulation to do precision tests of physics processes. 

This code was previously used with Geant4 v10.1.2 but not documented very well. We have migrated
to v10.7.4 and ROOT 6.24.2; both substantially advanced versions from what we started with, but
still not cutting edge (G4 is up to v11.1.1 and ROOT is up to v6.28.0 as of this writing). 


# Installing with G4 v10.7.4, ROOT 6.24.2, Ubuntu 20, gcc 9.4.0, standard=cxx14.


## Geant4 Install:

 - Download the source code: [v10.7.4](https://geant4.web.cern.ch/download/10.7.4.html)
 - make a build directory.
 - inside the build directory do:

```
cmake -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_INSTALL_DATA=ON -DCMAKE_CXX_STANDARD=14 -DCMAKE_PREFIX_PATH="/home/villaa/install/geant4/geant4-v10.7.4-install/lib/cmake/Geant4/" -DCMAKE_INSTALL_PREFIX=/home/villaa/install/geant4/geant4-v10.7.4-install /home/villaa/install/geant4/geant4-v10.7.4
```

Note the flags, if you miss the visualization, it will not be able to do visualization, if you
miss the install data flag it won't install the data (this happens at the make install step). If
you miss the forcing of the `cxx14` standard it will compile with a standard that might be
inconsistent with ROOT. The compilation will fail, hard. 

 - also inside the build directory do `$make -jN`, N is the number of parallel processors 
 - then do `$make install` (check that data directories are being downloaded!!)

 - make sure after you're installed to set the environmental variables like:

`source /home/villaa/install/geant4/geant4-v10.7.4-install/bin/geant4.sh` <br>
`export G4WORKDIR=/home/userrname/geant4bin`

The last one is to ensure that the `make install` process puts the binaries in a central location.
It's not strictly needed and could be a different directory.

## Installing ROOT:

 - If `cxx14` is not available for your operating system, you may have to try to install from source and force this version. For me I used an older version of ROOT precisely so this would not happen, then I could install by binary. 

 - go to here: [ROOT v6.24.2](https://root.cern/releases/release-62402/) get the Ubuntu 20 release `.tar.gz`. Unpack it and source the setup file like:

```
wget https://root.cern/download/root_v6.26.10.Linux-centos8-x86_64-gcc8.5.tar.gz
tar -xzvf root_v6.26.10.Linux-centos8-x86_64-gcc8.5.tar.gz
source root/bin/thisroot.sh # also available: thisroot.{csh,fish,bat}
```

 - check that it is compiled with `cxx14` by doing `root-config --features` and looking at the very first entry. 

## Compiling/installing NeutronReflectometry:

 - clone it from the git repo like: `git clone https://github.com/villano-lab/G4PhysicsTesting.git`. If there has not been a release past v0.9.1, go to the branch for this version of Geant: `feature/G4v10-7-4-upgrade`

 - make a "build" directory anywhere you like and run a cmake like:

```
cmake -DG4VIS_USE=ON -DCMAKE_CXX_STANDARD=14 -DGeant4_DIR=/home/villaa/install/geant4/geant4-v10.7.4/lib/Geant4-10.7.4/ /home/villaa/G4PhysicsTesting/
```

 - you can grep on visualization if you want to see if vis is on or off. It should be controlled by the `G4VIS_USE` supplied to `cmake` as above.

 - use the executable in that directory `NeutronReflectometry`

 - if the environmental variable `G4WORKDIR` has been set properly, it should also be available in `/home/username/geant4bin/bin/Linux-g++/`

`export G4WORKDIR=/home/userrname/geant4bin`
