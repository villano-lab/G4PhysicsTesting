/*==================NeutReflect_PhysicsList.hh============================= 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Class to define modular physics lists for use with Geant4.
              
======================================================================*/


#ifndef NeutReflect_PhysicsList_h
#define NeutReflect_PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"


class NeutReflect_PhysicsList : public G4VModularPhysicsList
{
public :

  NeutReflect_PhysicsList();
  ~NeutReflect_PhysicsList();

public :
  virtual void SetCuts();

protected :

private :

  G4int VerboseLevel;

};

#endif
