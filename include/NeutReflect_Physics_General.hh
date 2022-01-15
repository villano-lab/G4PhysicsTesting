/*==================NeutReflect_Physics_General.hh======================= 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Class to support general physics
               processes included in the .cc file.  The type inherits
	       from the G4VPhysicsConstructor class in Geant4.  Typically
	       physics lists are grouped by their similarity so that they
	       can be included modularly in an object which is a G4VModularPhysicsList.
	       This object is prototyped in the NeutReflect_PhysicsList.hh file of
	       this project code. 
              
======================================================================*/

#ifndef NeutReflect_Physics_General_h
#define NeutReflect_Physics_General_h 1

#include "globals.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4Decay.hh"


class NeutReflect_Physics_General : public G4VPhysicsConstructor
{
public:
  NeutReflect_Physics_General(const G4String& name = "general");
  virtual ~NeutReflect_Physics_General();

public:
  // This method will be invoked in the Construct() method.
  // each particle type will be instantiated
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();

protected:
  G4Decay fDecayProcess;
};

#endif
