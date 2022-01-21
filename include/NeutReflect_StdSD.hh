/*==================NeutReflect_StdSD.hh=================================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Class for specifying a sensitive detector in the NeutReflect 
               simulation.
              
======================================================================*/

#ifndef NeutReflect_StdSD_h
#define NeutReflect_StdSD_h 1

#include "G4VSensitiveDetector.hh"
#include "NeutReflect_StdHit.hh"

class G4Step;
class G4HCofThisEvent;
class NeutReflect_PrimaryGeneratorAction;

class NeutReflect_StdSD : public G4VSensitiveDetector
{
public:
  NeutReflect_StdSD(G4String, G4int panelNb);
  ~NeutReflect_StdSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  NeutReflect_StdHitsCollection* stdCollection;
  G4int genericNb;

  G4bool isNCap(G4Track* track,G4Step *aStep);
  NeutReflect_PrimaryGeneratorAction *Ngen;

};

#endif
