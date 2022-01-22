/*==================NeutReflect_StdSD.cc=================================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Code support for specifying a sensitive detector in the NeutReflect 
               simulation.  The variables, like positions and energies
	       are specified here for a given sensitive detector.
              
======================================================================*/

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "G4Neutron.hh"
#include "NeutReflect_StdSD.hh"
#include "NeutReflect_PrimaryGeneratorAction.hh"
extern G4int origevt;

NeutReflect_StdSD::NeutReflect_StdSD(G4String name, G4int genericNb):G4VSensitiveDetector(name), genericNb(genericNb)
{
  G4String HCname;// = "Collection";
  HCname = name;
  collectionName.insert(HCname);
  Ngen = (NeutReflect_PrimaryGeneratorAction*)(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
}
NeutReflect_StdSD::~NeutReflect_StdSD() {

}
void NeutReflect_StdSD::Initialize(G4HCofThisEvent* HCE)
{
  stdCollection = new NeutReflect_StdHitsCollection(SensitiveDetectorName,collectionName[0]);

  G4int HCID = -1;
  if(HCID < 0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    //    G4cout << ">>>>>> NeutReflect_StdSD::Initialize HCID = " << HCID << G4endl;
  }
  HCE->AddHitsCollection(HCID, stdCollection);
}
G4bool NeutReflect_StdSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
  G4double ke= aStep->GetPreStepPoint()->GetKineticEnergy();
  G4double pid;
  G4double edep=aStep->GetTotalEnergyDeposit();


  if(aStep->GetTrack()->GetDefinition()->GetParticleType() == "nucleus") {
    pid = (aStep->GetTrack()->GetDefinition()->GetPDGCharge()) + (1000 * aStep->GetTrack()->GetDefinition()->GetBaryonNumber());
  }
  else {
    pid = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  }

  NeutReflect_StdHit* newHit = new NeutReflect_StdHit();
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetDetNb(aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber());

  //G4cout << ">>>>>> ***** (NeutReflect_StdSD::ProcessHits) DetNb = " << aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber() << G4endl;
  newHit->SetEdep   (edep);
  newHit->SetPos    (aStep->GetPostStepPoint()->GetPosition());
  newHit->SetPName  (aStep->GetTrack()->GetDefinition()->GetParticleName());

  G4double* dataVector = newHit->GetData();

  dataVector[1-1]  = 0; // Is not set in this function
  dataVector[2-1]  = 1000 + genericNb; // Generic Detector Number
  dataVector[3-1]  = 1e5 * aStep->GetTrack()->GetTrackID() +
                           aStep->GetTrack()->GetCurrentStepNumber();
  dataVector[4-1]  = 1e5 * aStep->GetTrack()->GetParentID();
  dataVector[5-1]  = pid; // Particle ID
  dataVector[6-1]  = aStep->GetPreStepPoint()->GetKineticEnergy();
  dataVector[7-1]  = aStep->GetTotalEnergyDeposit();

  dataVector[8-1]  = aStep->GetPostStepPoint()->GetMomentum().x();
  dataVector[9-1]  = aStep->GetPostStepPoint()->GetMomentum().y();
  dataVector[10-1] = aStep->GetPostStepPoint()->GetMomentum().z();
  dataVector[11-1] = aStep->GetPostStepPoint()->GetPosition().x();
  dataVector[12-1] = aStep->GetPostStepPoint()->GetPosition().y();
  dataVector[13-1] = aStep->GetPostStepPoint()->GetPosition().z();
  dataVector[14-1] = aStep->GetPostStepPoint()->GetGlobalTime();

  dataVector[15-1] = aStep->GetPreStepPoint()->GetMomentum().x();
  dataVector[16-1] = aStep->GetPreStepPoint()->GetMomentum().y();
  dataVector[17-1] = aStep->GetPreStepPoint()->GetMomentum().z();
  dataVector[18-1] = aStep->GetPreStepPoint()->GetPosition().x();
  dataVector[19-1] = aStep->GetPreStepPoint()->GetPosition().y();
  dataVector[20-1] = aStep->GetPreStepPoint()->GetPosition().z();
  dataVector[21-1] = aStep->GetPreStepPoint()->GetGlobalTime();

  //get some process information
  G4Track* track = aStep->GetTrack();

  dataVector[22-1] = isNCap(track,aStep);

  /*
     1. EV (starts with 1)
     2. DT (detector #)
     3. TS (track and step: ttttsssss)
     4. PA (parent track  : tttt00000)
     5. TY (type. -ve =>  : -zzzaaa)
     6. E1 (KE of this track@1)
     7. D3 (energy deposition@3)
     8. X1 (position@1)
     9. Y1
    10. Z1
    11. X3 (position@3)
    12. Y3
    13. Z3 
    14. PX1 (momentum@1)
    15. PY1
    16. PZ1
    17. muxpos;
    18. muypos;
    19. muzpos;
    20. mutheta;
    21. muphi;
    22. muenergy;
 */

  stdCollection->insert(newHit);

  //  newHit->Print();
  newHit->Draw();

  return true;
}
void NeutReflect_StdSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  if (verboseLevel>0) {
    G4int NbHits = stdCollection->entries();     

    G4double TotEnergyDep = 0;

    for (G4int ii=0; ii<NbHits; ii++) {
      TotEnergyDep += (*stdCollection)[ii]->GetEdep();
    }

    G4cout << "--------> Total deposited energy in " 
	   << SensitiveDetectorName << " : " 
	   << std::setw(5) << G4BestUnit(TotEnergyDep,"Energy")
	   << "\n" << G4endl;

  }
}
G4bool NeutReflect_StdSD::isNCap(G4Track *track, G4Step *aStep){

  static const G4ParticleDefinition* theNeutron = G4Neutron::Definition();

  // Must be a neutron which has been stopped by a nucleus
  if (track->GetDefinition() != theNeutron) return false;
  if (track->GetTrackStatus() != fStopAndKill) return false;
  
  // Only tag if there were no neutron secondaries (identifies a capture)
  const G4TrackVector* secondaries = aStep->GetSecondary();
  size_t nsecs = secondaries->size();
  for(size_t i=0; i < nsecs; ++i) {    //quicker to search in reverse order
    if ((*secondaries)[nsecs-1-i]->GetDefinition() == theNeutron) return false;
  }

  return true;
}
