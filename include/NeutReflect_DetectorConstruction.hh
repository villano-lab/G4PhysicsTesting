/*==================NeutReflect_DetectorConstruction.hh=================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE:  Class prototype skeleton of geometry class for the NeutReflect
                simulation.
              
======================================================================*/

#ifndef NeutReflect_DetectorConstruction_H
#define NeutReflect_DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include <map>

class NeutReflect_DetectorConstructionMessenger;
class G4VPhysicalVolume;

// Sensitive detector classes can be defined here 
//class <sensitive detector>;

// -----------------Geometry Construction Class--------------------------

class NeutReflect_DetectorConstruction : public G4VUserDetectorConstruction
{

public:
  NeutReflect_DetectorConstruction(G4int design=-1, G4double diskR=12.7, G4double diskt=8.0, G4bool isBePure=false);
  ~NeutReflect_DetectorConstruction();

public:
  G4VPhysicalVolume* Construct(); //actual construction
                                  //of all components

  void UpdateGeometry();

  void SetConstructGenericGeometryBool(G4bool newVal)           {ConstructGenericGeometryBool = newVal;}
  void SetConstructGenericTrackerInt(G4int newVal) {ConstructGenericTrackerInt = newVal;}
  void SetConstructGenericSensitiveInt(G4int newVal) {ConstructGenericSensitiveInt = newVal;}

  G4bool GetConstructGenericGeometryBool()      {return ConstructGenericGeometryBool;}
  G4int GetConstructGenericTrackerInt() {return ConstructGenericTrackerInt;}
  G4int GetConstructGenericSensitiveInt() {return ConstructGenericSensitiveInt;}
  G4int GetDesignNo() {return DesignNo;}

  std::map<G4String,G4int> *GetSensitiveList() { return &NeutReflectCollName;}
  G4int    GetNSensitive() {return NeutReflectCollName.size();}
  G4String GetSensitive(G4int);

  //some gets
  G4bool GetBePure(){return UseBePure;}
  G4double GetBeR(){return BeR;}
  G4double GetBet(){return Bet;}
  G4double GetVesselR(){return VesselR;}
  G4double GetSphDetR(){return SphDetR;}
  G4double GetDiskDisplacement(){return DiskDisplacement;}

  void ChangeDetectorConfig(G4int design=0);
private:

  G4double world_x, world_y, world_z;
  G4VPhysicalVolume* physicalWorld;
  G4VPhysicalVolume* genericVol;

  //SD map
  std::map<G4String,G4int> NeutReflectCollName;

  //parameters
  G4double mfp_Pb_1_8; //mean free path in Lead of 1.8MeV gamma (units mm)

  //a useful rotation
  G4RotationMatrix *xrot;
  G4RotationMatrix *negxrot;
  G4RotationMatrix *nullrot;

  G4bool ConstructGenericGeometryBool;
  G4bool UseBePure;
  G4int ConstructGenericTrackerInt;
  G4int ConstructGenericSensitiveInt;
  G4int DesignNo;

  G4Material* exampleMat; //starts as vacuum
  G4Material* Pb;
  G4Material* Ge;
  G4Material* Si;
  G4Material* stainlessSteel;
  G4Material* BeO;
  G4Material* BePure;
  G4Material* Plastic;
  G4Material* Steel;
  G4Material* Air;
  G4Material* Argon;
  G4Material* Neon;
  G4Material* Aluminium;
  G4Material* shieldCuMat, *shieldPbMat;
  G4Material* stillHe,*MCHe;
  G4Material* helium;
  G4Material* GasHe,*Gas3He;
  G4Material* Gas3He6atm,*Gas3He10atm;

  G4double BeR;
  G4double Bet;
  G4double VesselR;
  G4double SphDetR;
  G4double DiskDisplacement;

  NeutReflect_DetectorConstructionMessenger* detMessenger; 

  G4Region* DetectorRegion;

  void DefineMaterials();
  void ConstructTracking();
  void ConstructGenericSensitive(G4LogicalVolume*,G4String);
  void ConstructNULL(G4VPhysicalVolume*);
  void ConstructDiag0(G4VPhysicalVolume*);
  void ConstructDiag1(G4VPhysicalVolume*);
  void ConstructDiag2(G4VPhysicalVolume*);
  void ConstructDiag3(G4VPhysicalVolume*);
  void ConstructDiag4(G4VPhysicalVolume*);
  void ConstructDesign0(G4VPhysicalVolume*,G4double tl=6.0);
  void ConstructDesign1(G4VPhysicalVolume*,G4double tl=6.0,
					   G4double tu=100.0,
					   G4double dr=30.0,
					   G4double g1=10.0,
					   G4double g2=10.0,
					   G4double hr=200.0,
					   G4double rr=30.0);
  void ConstructDesign2(G4VPhysicalVolume*,G4double tl=6.0,
                                           G4double gedist=400.0,
                                           G4double gerad=76.2,
 					   G4double gethk=25.4+12.7,
                                           G4double beorad=20.0);
  void ConstructDesign3(G4VPhysicalVolume*,G4double tl=6.0,
					   G4double dpipe=2.5,
					   G4double tseteelplate=40.0,
					   G4double tsrcholder=6.35,
					   G4double lBeO=400.0,
					   G4double ldet=186.2,
					   G4double ddet=152.4);
  void ConstructSourcePlug(G4VPhysicalVolume*,G4double xpos=0.0,
					      G4double ypos=0.0,
 					      G4double zpos=0.0);
  void ConstructStainlessPlug(G4VPhysicalVolume*,G4double xpos=0.0,
					      G4double ypos=0.0,
 					      G4double zpos=0.0);
  void ConstructBePlug(G4VPhysicalVolume*,G4double xpos=0.0,
					  G4double ypos=0.0,
					  G4double zpos=-4.0-(6.35/2.0));
  void ConstructGeDet(G4VPhysicalVolume*,G4double dist=400.0);
  void ConstructGeDetSoudan(G4VPhysicalVolume*,G4double dist=400.0);
  void ConstructSiDet(G4VPhysicalVolume*,G4double dist=400.0);
  void ConstructR68SiDet(G4VPhysicalVolume*,G4double dist=400.0);
  void ConstructSNOLABDet(G4VPhysicalVolume*,G4String name, G4String mat="Ge",G4double dist=400.0);
  void ConstructArDet(G4VPhysicalVolume*,G4double dist=400.0);
  void ConstructNeDet(G4VPhysicalVolume*,G4double dist=400.0);
  void ConstructGenericDet(G4VPhysicalVolume*,G4Material*,G4String name="gendet",G4double dist=400.0);
  void ConstructSimpleVessel(G4VPhysicalVolume*,G4double vthk=400.0,
						G4double vrad=146.2,
						G4double vh=2000.0);
  void ConstructPhotoNSourceBox(G4VPhysicalVolume*,G4double xpos=0.0,
					      G4double ypos=0.0,
 					      G4double zpos=0.0,
					      G4String file8020="8020_medres.txt",
                                              G4int bshelf=0,
                                              G4int gshelf=1,
					      G4bool dolead=true,
					      G4bool reston=true);
  void ConstructPhotoNSourceBoxNested(G4VPhysicalVolume*,G4double xpos=0.0,
					      G4double ypos=0.0,
 					      G4double zpos=0.0,
                                              G4int bshelf=0,
                                              G4int gshelf=1);
  void Construct3HeTube(G4VPhysicalVolume*);
  void ConstructL3HeDesign1(G4VPhysicalVolume*);
  

};

#endif
