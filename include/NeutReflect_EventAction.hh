/*==================NeutReflect_EventAction.hh============================= 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE:  Class to support the track actions carried out
                when an "event" occurs in the simulation.
              
======================================================================*/

#ifndef NeutReflect_EventAction_h
#define NeutReflect_EventAction_h 1

#include "G4UserEventAction.hh"

class NeutReflect_RunAction;
class NeutReflect_DataStorage;
class NeutReflect_DetectorConstruction;

// ------------------------------------------

class NeutReflect_EventAction : public G4UserEventAction
{
public :

  NeutReflect_EventAction();
  NeutReflect_EventAction(NeutReflect_RunAction*);
  ~NeutReflect_EventAction();

public :

  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

private :


  G4bool drawEvent, saveEvent;

  //data output class. not used for initial encarnation
  NeutReflect_DataStorage* dataOut;
  NeutReflect_DetectorConstruction* NeutReflect_Detector;
  NeutReflect_RunAction* pRunAction;

};

#endif
