/*==================NeutReflect_DetectorConstructionMessenger.hh============ 
   
      PROGRAMMER:  Anthony Villano 10/04/11

      UPDATES:      
       

      PURPOSE: Class to supply an interactive interface to the DetectorConstruction
      	       class and make the setting of the geometrical parameters
	       able to be done in an interactive session or script. 
              
              
======================================================================*/
#ifndef NeutReflect_DetectorConstructionMessenger_h
#define NeutReflect_DetectorConstructionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

// ------------------------------------------------

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

class NeutReflect_DetectorConstruction;

// ------------------------------------------------

class NeutReflect_DetectorConstructionMessenger: public G4UImessenger
{
public:
  NeutReflect_DetectorConstructionMessenger(NeutReflect_DetectorConstruction* );
  ~NeutReflect_DetectorConstructionMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  NeutReflect_DetectorConstruction*       pDetCons;
  G4UIdirectory*        NeutReflect_DetConfigChange;

  // Run Element Activation

  G4UIcmdWithAnInteger*       setNewDetConfig;
  G4UIcmdWithoutParameter*       updateGeom;
};
#endif
