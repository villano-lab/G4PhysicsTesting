/*==================NeutReflect_DetectorConstructionMessenger.cc=========== 
   
      PROGRAMMER:  Anthony Villano 10/20/12

      UPDATES:      
       

      PURPOSE: Code supporting the DetectorConstructionMessenger class including
               specific routines which run DetectorConstruction methods with various
	       input depending on the values passed via Geant4 macro.
              
              
======================================================================*/
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"


#include "G4RunManager.hh"

#include "NeutReflect_DetectorConstructionMessenger.hh"
#include "NeutReflect_DetectorConstruction.hh"
#include <iostream>
#include <sstream>
#include <string>

NeutReflect_DetectorConstructionMessenger::NeutReflect_DetectorConstructionMessenger(NeutReflect_DetectorConstruction* NRgeom):pDetCons(NRgeom)
{ 
  // Create the run directory
  NeutReflect_DetConfigChange = new G4UIdirectory("/det/NeutReflect/");
  NeutReflect_DetConfigChange->SetGuidance("NeutReflect specific detector controls.");


  // Set new configuration 
  setNewDetConfig = new G4UIcmdWithAnInteger("/det/NeutReflect/Configuration",this);
  setNewDetConfig->SetGuidance("Set the detector configuration with integer code.");
  setNewDetConfig->SetParameterName("fname",true);
  setNewDetConfig->SetDefaultValue(0);
  setNewDetConfig->AvailableForStates(G4State_PreInit,G4State_Idle);

  // Pure geometry update 
  updateGeom = new G4UIcmdWithoutParameter("/det/NeutReflect/update",this);
  updateGeom->SetGuidance("Rebuild all geometry.");
  updateGeom->SetGuidance("This command MUST be applied before \"beamOn\" ");
  updateGeom->SetGuidance("if you changed geometrical value(s).");
  updateGeom->AvailableForStates(G4State_Idle);


}
NeutReflect_DetectorConstructionMessenger::~NeutReflect_DetectorConstructionMessenger()
{

}
void NeutReflect_DetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 

  if( command == setNewDetConfig ){
    std::istringstream valstream(newValue);
    G4int intval;
    valstream >> intval;
    pDetCons->ChangeDetectorConfig(intval);
    pDetCons->UpdateGeometry();
  }
  if( command == updateGeom ){
    pDetCons->UpdateGeometry();
  }
 
}
