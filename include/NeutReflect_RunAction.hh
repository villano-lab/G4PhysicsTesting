/*==================NeutReflect_RunAction.hh================================= 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Class to keep track of runs as instantiated from a Geant4
               macro.  Every /run/beamOn n command begins a new run of
	       n thrown events.
              
======================================================================*/

#ifndef NeutReflect_RunAction_h
#define NeutReflect_RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class NeutReflect_DataStorage;
class NeutReflect_RunActionMessenger;
class NeutReflect_PrimaryGeneratorAction;

class NeutReflect_RunAction : public G4UserRunAction
{
public :

  NeutReflect_RunAction(NeutReflect_PrimaryGeneratorAction*);
  ~NeutReflect_RunAction();

public :

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

  void SetDataFileNamePrefix(G4String nPrefix) {DataFileNamePrefix = nPrefix;}
  void SetUseRandomTag(G4bool usetag) {UseRandomTag = usetag;}
  void SetOutputDataToFile(G4bool dataout)     {OutputDataToFile = dataout;}
  void SetDrawEventCmd(G4bool drawBool)        {DrawEventCmd = drawBool;}
  void SetSaltPillOutCmd(G4bool value)         {SaltPillOutCmd = value;}
  void SetAutoSeed (const G4bool value)        {autoSeed    = value;}

  G4String GetDataFileNamePrefix() {return DataFileNamePrefix;}
  G4bool GetUseRandomTag() {return UseRandomTag;}
  G4bool GetOutputDataToFile()   {return OutputDataToFile;}
  G4bool GetDrawEventCmd()       {return DrawEventCmd;}
  G4bool GetSaltPillOutCmd()     {return SaltPillOutCmd;}
  G4int  GetRunNumber()          {return runN;}

  NeutReflect_DataStorage* dataOut;

private :

  NeutReflect_RunActionMessenger* runMessenger; 

  G4bool   autoSeed;
  G4String DataFileNamePrefix;
  G4bool   UseRandomTag;
  G4bool   DrawEventCmd;
  G4bool   SaltPillOutCmd;
  G4bool   OutputDataToFile;
  G4int    runN;
  long     randSeed;
  NeutReflect_PrimaryGeneratorAction* thisgenerator;

};

// ------------------------------------------


#endif
