/*==================NeutReflect_PrimaryGeneratorAction.hh================= 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE: Class which defines the generators for the NeutReflect
               simulation.
              
======================================================================*/



#ifndef NeutReflect_PrimaryGeneratorAction_h
#define NeutReflect_PrimaryGeneratorAction_h 1

#include "G4ParticleTable.hh"
#include "globals.hh"

#include "G4VUserPrimaryGeneratorAction.hh"

#include <vector>

//root stuff
#include "TFile.h"
#include "TTree.h"
#include "Rtypes.h"

// ------------------------------------------------

class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class NeutReflect_DetectorConstruction;

// ------------------------------------------------
class NeutReflect_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

public:

  NeutReflect_PrimaryGeneratorAction(G4double bias=1.0, G4String dec="88Y", G4String rootfile="NULL");
  ~NeutReflect_PrimaryGeneratorAction();

   void SetSourceType(G4String type){ sourceType=type;}

   G4String GetSourceType() {return sourceType;}
   G4double GetFinalPathLength() {return finalpathlength;}
 
public:

  void GeneratePrimaries(G4Event* anEvent);

private:

  NeutReflect_DetectorConstruction *Ndet; 

  //root file reading
  TFile *f;
  TTree *t;
  int rootnevent;  //keep track of which event we're on
  Long64_t ncascade;
  double time[1000],Eg[1000];

  //kinematics
  G4double mbe,mn,mnbar,melec,deln,delbe; 

  //cross section biasing -- in cm sized samples it's 
  //probably safe to bias as much as 1000 times without
  //worrying about attenuation or correction to the gamma ratio
  G4double xnbias;

  G4String sourceType;
  G4String decayer;
  G4String rootfile;
  G4double R,thk,d;
  G4ParticleGun* particleGun;
  G4RotationMatrix *xrot;
  G4double finalpathlength;
  //G4GeneralParticleSource* particleGun;

  //helpful functions
  std::vector<G4double> GenerateRandomDirection();
  std::vector<G4double> GenerateGammaEnergy(G4String spectrum="88Y");
  void GenerateGammaLevel88Sr(std::vector<G4double> &gammas,G4int lev=0);
  void GenerateGammaLevel124Te(std::vector<G4double> &gammas,G4int lev=0);
  std::vector<G4double> GenerateDiskPathlength(std::vector<G4double> angles);
  std::vector<G4double> GenerateShellPathlength(std::vector<G4double> angles);
  G4double GenerateGammaNEnergy(G4double costhet,G4double egamma,G4double mrec, G4double Q);
  G4bool DidGammaNHappen(G4double pathlength,G4double energy);


};

#endif
