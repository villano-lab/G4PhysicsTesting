/*==================NeutReflect_RunActionMessenger.hh===================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Class to supply an interactive interface to the RunAction
               class.  Examples of tasks include setting the file name
	       prefix for the output. This class will supply commands
	       accessible in standard Geant4 Macros
      
              
======================================================================*/

#ifndef NeutReflect_RunActionMessenger_h
#define NeutReflect_RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

// ------------------------------------------------

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

class NeutReflect_RunAction;

// ------------------------------------------------

class NeutReflect_RunActionMessenger: public G4UImessenger
{
public:
  NeutReflect_RunActionMessenger(NeutReflect_RunAction* );
  ~NeutReflect_RunActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  NeutReflect_RunAction*       pRunAction;
  G4UIdirectory*        NeutReflect_RunDir;

  // Run Element Activation

  G4UIcmdWithABool*       setDrawEventCmd;
  G4UIcmdWithABool*       setUseRandomTagCmd;
  G4UIcmdWithAString*       setRunFileName;
};

#endif
