* *be sure to have the Geant4 environment sourced* this is the step I always forget, for me it was something like:

. ~/install/geant4/geant4.10.07.p03-install/bin/geant4.sh    OR

source ~/install/geant4/geant4.10.07.p03-install/bin/geant4.sh


* Make a separate directory and call it G4PhysicsTest_build or something

* run cmake like:

cmake -DGeant4_DIR=/Users/villaa/install/geant4/geant4.10.07.p03-install -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_QT=ON /Users/villaa/G4PhysicsTesting

The first variable is the *full path* to your Geant4 installation; the graphics variables probably don't work--issue for another time; and the last thing is the **full path** to the G4PhysicsTesting repository. 

* type `make`

* if it succeeds you're either super awesome, lucky, or both. 
