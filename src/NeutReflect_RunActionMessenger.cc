/*==================NeutReflect_RunActionMessenger.cc====================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE: Code supporting the RunActionMessenger class including
               specific routines which run RunAction methods with various
	       input depending on the values passed via Geant4 macro.
              
======================================================================*/

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

#include "G4RunManager.hh"

#include "NeutReflect_RunActionMessenger.hh"
#include "NeutReflect_RunAction.hh"


NeutReflect_RunActionMessenger::NeutReflect_RunActionMessenger(NeutReflect_RunAction* NeutReflect_Run):pRunAction(NeutReflect_Run)
{ 
  // Create the run directory
  NeutReflect_RunDir = new G4UIdirectory("/run/NeutReflect/");
  NeutReflect_RunDir->SetGuidance("NeutReflect specific run controls.");

  //  run directory already exists

  // Set autoSeed command
  //setAutoSeedCmd = new G4UIcmdWithABool("/run/autoSeed",this);
  //setAutoSeedCmd->SetGuidance("Switch on/off time-based random seeds");
  //setAutoSeedCmd->SetGuidance(" true: run seeds determined by system time");
  //setAutoSeedCmd->SetGuidance("false: use command 'random/resetEngineFrom'");
  //setAutoSeedCmd->SetGuidance("Default = false");
  //setAutoSeedCmd->SetParameterName("autoSeed", false);
  //setAutoSeedCmd->AvailableForStates(G4State_Idle);

  // Set OutputDataToFile
  //setOutputDataToFile = new G4UIcmdWithABool("/run/writedata",this);
  //setOutputDataToFile->SetGuidance("Output results to file.");
  //setOutputDataToFile->SetParameterName("dataout",true,true); 
  // void SetParameterName(const char * theName,G4bool omittable,
  //                      G4bool currentAsDefault=false);
  //setOutputDataToFile->SetDefaultValue(true);
  //setOutputDataToFile->AvailableForStates(G4State_PreInit,G4State_Idle);

  // Set file name
  setRunFileName = new G4UIcmdWithAString("/run/NeutReflect/OFPrefix",this);
  setRunFileName->SetGuidance("Set the name of the output files.");
  setRunFileName->SetParameterName("fname",true);
  setRunFileName->SetDefaultValue("Sim");
  setRunFileName->AvailableForStates(G4State_PreInit,G4State_Idle);

  // Set drawEvent
  //  the  G4UIdirectory("/NeutReflect/"); directory is defined in DetectorMessenger
  setDrawEventCmd = new G4UIcmdWithABool("/run/NeutReflect/Draw",this);
  setDrawEventCmd->SetGuidance("Set drawFlag to Draw an event. (default = true)");
  setDrawEventCmd->SetParameterName("Draw", false);
  setDrawEventCmd->SetDefaultValue(true);

  // set or unset the UserRandomTag (Default true)
  //  the  G4UIdirectory("/NeutReflect/"); directory is defined in DetectorMessenger
  setUseRandomTagCmd = new G4UIcmdWithABool("/run/NeutReflect/UseRandomTag",this);
  setUseRandomTagCmd->SetGuidance("Set UseRandomTag for filenames. (default = true)");
  setUseRandomTagCmd->SetParameterName("UseRandomTag", false);
  setUseRandomTagCmd->SetDefaultValue(true);
}
NeutReflect_RunActionMessenger::~NeutReflect_RunActionMessenger()
{
   delete setDrawEventCmd;     
}
void NeutReflect_RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 


  if( command == setDrawEventCmd ) {
    if(setDrawEventCmd->GetNewBoolValue(newValue)) {
      G4cout << "Turning ON event drawing for all events." << G4endl; }
    else {
      G4cout << "Turning OFF event drawing for all events." << G4endl; }
    pRunAction->SetDrawEventCmd(setDrawEventCmd->GetNewBoolValue(newValue));
  }

  if( command == setUseRandomTagCmd ) {
    if(setUseRandomTagCmd->GetNewBoolValue(newValue)) {
      G4cout << "Turning ON random filename tags." << G4endl; }
    else {
      G4cout << "Turning OFF random filename tags." << G4endl; }
    pRunAction->SetUseRandomTag(setUseRandomTagCmd->GetNewBoolValue(newValue));
  }

  if( command == setRunFileName ){
    pRunAction->SetDataFileNamePrefix(newValue);
  }

  
}
