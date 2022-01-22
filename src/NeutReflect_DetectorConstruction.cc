/*==================NeutReflect_DetectorConstruction.cc==================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Code support for detector construction in the NeutReflect
               simulation.  It is attempted to make geometry defined
	       somewhat generally to accomodate design changes after
	       this files inception date. 
              
======================================================================*/

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4ExtrudedSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "NeutReflect_StdSD.hh"
#include "NeutReflect_DetectorConstructionMessenger.hh"
#include "NeutReflect_DetectorConstruction.hh"

//--------------Geometry Constructor ---------------------

NeutReflect_DetectorConstruction::NeutReflect_DetectorConstruction(G4int design, G4double diskR, G4double diskt, G4bool isBePure):
UseBePure(isBePure), //Beoxide or Be pure default BeO
BeR(diskR),
Bet(diskt),
VesselR(146.2), //I hate that this is basically hard-coded
SphDetR(2*m), 
DiskDisplacement(-3.175)
{
	
	// -------- The World ---------
	//world_x = 16.01*m; world_y = 16.01*m; world_z = 9.0*m;
	world_x = 3.0*m; world_y = 3.0*m; world_z = 3.0*m;

        //set some parameters
        mfp_Pb_1_8=18.6;

        //a useful rotation
	G4ThreeVector row1 = G4ThreeVector(1.0,0,0);
	G4ThreeVector row2 = G4ThreeVector(0,0,1.0);
	G4ThreeVector row3 = G4ThreeVector(0,-1.0,0);
	xrot = new G4RotationMatrix(row1,row2,row3);
	xrot->setRows(row1,row2,row3);

        //invers of a useful rotation
	row1 = G4ThreeVector(1.0,0,0);
	row2 = G4ThreeVector(0,0,-1.0);
	row3 = G4ThreeVector(0,+1.0,0);
	negxrot = new G4RotationMatrix(row1,row2,row3);
	negxrot->setRows(row1,row2,row3);

        //a null rotation
	row1 = G4ThreeVector(1.0,0,0);
	row2 = G4ThreeVector(0,1.0,0);
	row3 = G4ThreeVector(0,0,1.0);
	nullrot = new G4RotationMatrix(row1,row2,row3);
	nullrot->setRows(row1,row2,row3);
	
	
	// Create commands for interactive definition of the detector
	//detectorMessenger = new NeutReflect_DetectorConstructionMessenger(this);
	
	
	//turn on and off different aspects
	ConstructGenericGeometryBool = true;
	//tracking aspects
	ConstructGenericTrackerInt=0;
        //sensitive detectors
        ConstructGenericSensitiveInt=1;
	
	
	// ---------Material Definition--------------
	DefineMaterials();

        //set the design number	
        DesignNo=design;

	// ---------Detector Names--------------
	//NeutReflectCollName = new char*[30];
	//NeutReflectCollName[0]  = "NeutReflect_Detector1";
	
        //make messenger 
        detMessenger = new NeutReflect_DetectorConstructionMessenger(this);  
}

// ---------------Geometry Destructor------------------

NeutReflect_DetectorConstruction::~NeutReflect_DetectorConstruction()
{
	// Don't have anything here to delete
}

// ---------------Materials Definition------------------------

void NeutReflect_DetectorConstruction::DefineMaterials()
{
	
	G4String name, symbol;    
	G4double a, z,iz, density;
	
	G4int nel,ncomponents, natoms;
	//G4double fractionmass;
	G4double temperature, pressure;
	
	
        //Pb
        a = 207.2*g/mole;
        G4Element* elPb = new G4Element(name="Lead",symbol="Pb",z=82.,a);
        density = 11.35*g/cm3;
        Pb = new G4Material(name="Pb",density,ncomponents=1);
        Pb->AddElement(elPb,1);      
 
 
        //Copper
        a=63.5460*g/mole;
        density=8.920*g/cm3;
        G4Element* elCu = new G4Element(name="Copper", symbol="Cu", z=29., a);
        G4Material* Cu = new G4Material(name="Cu", density, ncomponents=1);
        Cu->AddElement(elCu, 1); 


        //Si
        a = 28.09*g/mole;
        density = 2.330*g/cm3;
        G4Element* elSi = new G4Element(name="Silicon",symbol="Si",z=14.,a);
        Si=new G4Material(name="Silicon",density, ncomponents=1);
        Si->AddElement(elSi,1);

        //Ge
        a=72.61*g/mole;
        density=5.323*g/cm3;
        G4Element* elGe = new G4Element
           (name="Germanium", symbol="Ge", z=32., a);// doesn't work with inel
        Ge = new G4Material(name="Ge", density, ncomponents=1);
        Ge->AddElement(elGe, 1);

        //Al
        a=26.9815*g/mole;
        G4Element* elAl=new G4Element(name="Aluminium",symbol="Al",z=13.,a);

        //C
        a = 12.011*g/mole;
        G4Element* elC=new G4Element(name="Carbon",symbol="C",z=6.,a);
          
        //H
        a=1.01*g/mole;
        G4Element* elH=new G4Element(name="Hydrogen",symbol="H",z=1.,a);

        //Poly
        density = 0.935*g/cm3;
        G4Material* poly=new G4Material(name="Poly", density,ncomponents=2);
        poly->AddElement(elH,natoms=2);
        poly->AddElement(elC,natoms=1);

        //Scintillator
        density=1.032*g/cm3;
        G4Material* Scint=new G4Material(name="Scintillator",density,ncomponents=2);
        Scint->AddElement(elH, natoms=11);
        Scint->AddElement(elC, natoms=10);

        // Aluminium
        Aluminium = new G4Material(name="Aluminium", density = 2.700*g/cm3, ncomponents=1);
        Aluminium->AddElement(elAl, natoms=1);

        //Air
        a = 14.01*g/mole;
        G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
        a = 16.00*g/mole;
        G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);
        density = 1.29*mg/cm3;
        //STP is default so don't have to supply temp and press
        Air = new G4Material(name="Air", density, nel=2, kStateGas);
        Air->AddElement(elN, 0.7);
        Air->AddElement(elO, 0.3);

        //Argon
        a = 39.948*g/mole;
        G4Element* elAr = new G4Element(name="Argon", symbol="Ar", iz=18.0, a);
        density = 1.784*mg/cm3;
        //STP is default so don't have to supply temp and press
        Argon = new G4Material(name="Argon", density, nel=1, kStateGas);
        Argon->AddElement(elAr, 1.0);
 
        //Neon
        a = 20.18*g/mole;
        G4Element* elNe = new G4Element(name="Neon", symbol="Ne", iz=10.0, a);
        density = 0.89994*mg/cm3;
        //STP is default so don't have to supply temp and press
        Neon = new G4Material(name="Neon", density, nel=1, kStateGas);
        Neon->AddElement(elNe, 1.0);

        //Steel (i.e. iron)
        a=55.845*g/mole;
        density=7.874*g/cm3;
        G4Element* elFe=new G4Element(name="Iron", symbol="Fe", z=26., a);
        Steel=new G4Material(name="Steel",density,ncomponents=1);
        Steel->AddElement(elFe, 1);

        //------------------------------------------------
        //Define new Elements for BeO
        //-----------------------------------------------
	G4double aBe,aO;
	
	//Be
	aBe= 9.012182*g/mole;
	G4Element* elBe=new G4Element(name="Beryllium",symbol="Be",z=4.,aBe);

	//O
	aO= 15.9994*g/mole;
	//G4Element* elO=new G4Element(name="Oxygen",symbol="O",z=8.,aO);

        //molar mass of Be0 is 25 g/mol and cross section is about 1.5 mb for gammas in
        //the right energy range, together this means that the mfp for this
        //interaction is about 184 meters. So the correctiona is about 10^-3 % for the
        //next term of the attenuating exponential. Increasing the cross section by a factor
        //of 1000 induces an approx 5% error in this next term -- this is as far as I would
        //bias it for cm sized samples
        density=3.02*g/cm3;
        BeO = new G4Material(name="BeO",density,ncomponents=2); 
        BeO->AddElement(elBe, natoms=1);
        BeO->AddElement(elO, natoms=1);

        //also interested in pure Be
        density=1.85*g/cm3;
        BePure = new G4Material(name="BePure",density,ncomponents=1); 
        BePure->AddElement(elBe, natoms=1);
	// ------------------------------------------------
	//Define new Elements for Stainless Steel Type 304
	// ------------------------------------------------
	G4double aC,aFe,aCr,aMn,aNi,aP,aSi,aS;

	//C
	aC = 12.011*g/mole;
	//G4Element* elC=new G4Element(name="Carbon",symbol="C",z=6.,aC);
	
        //Fe
	aFe=55.845*g/mole;
	//G4Element* elFe=new G4Element(name="Iron",symbol="Fe",z=26.,aFe);

	//Cr
	aCr=51.996*g/mole;
	G4Element* elCr=new G4Element(name="Chromium",symbol="Cr",z=24.,aCr);

	//Ni
	aNi=58.693*g/mole;
	G4Element* elNi=new G4Element(name="Nickel",symbol="Ni",z=28.,aNi);

	//S
	aS=32.065*g/mole;
	G4Element* elS=new G4Element(name="Sulfur",symbol="S",z=16.,aS);

	//Si
	aSi=28.086*g/mole;
	//G4Element* elSi=new G4Element(name="Silicon",symbol="Si",z=14.,aSi);

	//P
	aP=30.974*g/mole;
	G4Element* elP=new G4Element(name="Phosphorus",symbol="P",z=15.,aP);

	//Mn
	aMn=54.938*g/mole;
	G4Element* elMn=new G4Element(name="Manganese",symbol="Mn",z=25.,aMn);
	// Vacuum
	density     = universe_mean_density;   //from PhysicalConstants.h
	                                       //which is loaded in some precompiled
					       // library, CLEP maybe?
	pressure    = 1.0e-19*pascal; 
	temperature = 2.73*kelvin;
	G4Material* Vacuum = new G4Material("Vacuum", z=1.0, a=1.01*g/mole, density,
										kStateGas, temperature, pressure);

	// Stainless Steel Type 304 Alloy 
	G4double percFe,percC=0.04,percCr=19.0,percMn=1.0,percNi=9.25,percP=0.0225,percSi=0.5,percS=0.015;
	G4double massfracFe,massfracC,massfracCr,massfracMn,massfracNi,massfracP,massfracSi,massfracS,masstot;
	percFe=100.0-(percC+percCr+percMn+percNi+percP+percSi+percS);
	masstot=(percFe*aFe+percC*aC+percCr*aCr+percMn*aMn+percNi*aNi+percP*aP+percSi*aSi+percS*aS);
	massfracFe=percFe*aFe/masstot;
	massfracC=percC*aC/masstot;
	massfracCr=percCr*aCr/masstot;
	massfracMn=percMn*aMn/masstot;
	massfracNi=percNi*aNi/masstot;
	massfracP=percP*aP/masstot;
	massfracSi=percSi*aSi/masstot;
	massfracS=percS*aS/masstot;
	G4double lbtog=453.59237,invin3toinvcm3=1/pow(2.54,3.0);
	stainlessSteel = new G4Material("Stainless", density = 0.285*lbtog*invin3toinvcm3*g/cm3 , ncomponents=8);
	stainlessSteel->AddElement(elFe,massfracFe);
	stainlessSteel->AddElement(elC,massfracC);
	stainlessSteel->AddElement(elCr,massfracCr);
	stainlessSteel->AddElement(elMn,massfracMn);
	stainlessSteel->AddElement(elNi,massfracNi);
	stainlessSteel->AddElement(elP,massfracP);
	stainlessSteel->AddElement(elSi,massfracSi);
	stainlessSteel->AddElement(elS,massfracS);
	
	// ------------------------------------------------
	//Define new Elements for "high strength plastic" 
	// ------------------------------------------------
	//G4double aTEST; 
        Plastic = poly;


	// Assign materials to some of the major components

        exampleMat = Vacuum; //change to air at STP
	
	// ------------------------------------------------
	
	
}
// --------------Construct NeutReflect Setup-------------------------------

void NeutReflect_DetectorConstruction::ChangeDetectorConfig(G4int design)
{
        G4GeometryManager::GetInstance()->OpenGeometry();
        G4PhysicalVolumeStore::GetInstance()->Clean();
        G4LogicalVolumeStore::GetInstance()->Clean();
        G4SolidStore::GetInstance()->Clean();

	NeutReflectCollName.clear();
        DesignNo=design;
	/*G4SDManager* SDman = G4SDManager::GetSDMpointer();

        G4String thickSDname = "GFcompare_Detector1";
        G4String thinSDname = "GFcompare_Detector2";
        GFcompare_StdSD* SD1 = (GFcompare_StdSD*) SDman->FindSensitiveDetector(thickSDname);
        GFcompare_StdSD* SD2 = (GFcompare_StdSD*) SDman->FindSensitiveDetector(thinSDname);
	*/
	/*if(SD1){
	  delete SD1;
	}
	if(SD2){
	  delete SD2;
	}*/

        return;
}
G4VPhysicalVolume* NeutReflect_DetectorConstruction::Construct()
{
	// ------------ Construct the Physical world ---------------
	
	G4Box* solidWorld = new G4Box("world_S", world_x, world_y, world_z);
	G4LogicalVolume*  logicalWorld = new G4LogicalVolume(solidWorld,  // The solid
									 exampleMat, // Material
									 "world_L",  // Name
									 0,0,0);
	// Physical volume
	physicalWorld = new G4PVPlacement(0,// no rotation
					  G4ThreeVector(), // at (0, 0, 0)
					  "world_P",       // name (using the second constructor)
					  logicalWorld,    // the logical volume to use
					  NULL,            // the mother volume
					  false,           // no boolean operation
					  0);              // copy number

	// Visualization attributes
	G4VisAttributes* VisAttWorld = new G4VisAttributes(G4Colour(0.8,1.0,1.0));
	VisAttWorld->SetForceWireframe(true);  //I want a Wireframe of the total volume
	logicalWorld->SetVisAttributes(VisAttWorld);
	logicalWorld->SetVisAttributes(G4VisAttributes::Invisible);
	
	// --------- End Construct the Physical world --------------
	
	//------------------------------------------------ 
	// Construct other detector components 
	//------------------------------------------------ 
        //Parameters
        G4double dist=400.0;
        G4double tl=6.0,tu=100.0,dr=30.0,g1=10.0,g2=10.0,hr=200.0,rr=30.0;
        G4double dpipe=5.0,tsteelplate=20.0,tsrcholder=6.35,lBeO=400.0,ldet=186.2,ddet=152.4;
        //G4double dpipe=5.0,tsteelplate=6.35,tsrcholder=6.35,lBeO=400.0,ldet=186.2,ddet=152.4;
        //start by constructing Ge detector
        if(DesignNo>-2){ //use below -2 for diagnostic geometries
          ConstructGeDet(physicalWorld,dist);
          ConstructSimpleVessel(physicalWorld,40.0,76.2+dr+g2+rr);
        }

        if(DesignNo==0){
          //construct Design0
          ConstructDesign0(physicalWorld,tl);
        }
        else if(DesignNo==1){
          //construct Design1
          ConstructDesign1(physicalWorld,tl,tu,dr,g1,g2,hr,rr);
        }
        else if(DesignNo==2){
          //construct Design2
          ConstructDesign2(physicalWorld,tl,dist);
        }
        else if(DesignNo==3){
          //construct Design3
          ConstructDesign3(physicalWorld,tl,dpipe,tsteelplate,tsrcholder,lBeO,ldet,ddet);
        }
        else if(DesignNo==-1){
          //construct NULL
          ConstructNULL(physicalWorld);
        }
        else if(DesignNo==-2){
          //construct Diag0 
          ConstructDiag0(physicalWorld);
        }
        else if(DesignNo==-3){
          //construct Diag1 
          ConstructDiag1(physicalWorld);
        }
        else if(DesignNo==-4){
          //construct Diag2 
          ConstructDiag2(physicalWorld);
        }
        else if(DesignNo==-5){
          //construct Diag3
	  //same as Diag2 but with no lead stack
          ConstructDiag3(physicalWorld);
        }
        else if(DesignNo==-6){
          //construct Diag4
	  //same as Diag3 but with nothing but Be disk
          ConstructDiag4(physicalWorld);
        }
        else if(DesignNo==-7){
          //construct Diag5
	  //just a Ge detector centered at zero
          ConstructGeDet(physicalWorld,0.0);
        }
        else if(DesignNo==-8){
          //construct Diag6
	  //just a Si detector centered at zero
          ConstructSiDet(physicalWorld,0.0);
        }
        else if(DesignNo==-9){
          //construct Diag7
	  //just a Ar detector centered at zero
          ConstructArDet(physicalWorld,0.0);
        }
        else if(DesignNo==-10){
          //construct Diag8
	  //just a Ge/Si tower of SNOLAB dets 
	  //for cascade generation in the Ge detectors ONLY
	  G4double spacing = 33.0+10.16; //use 0.4in spacing
	  //spacing negative switches order of towers
	  spacing *= -1;
          ConstructSNOLABDet(physicalWorld,"Ge1","Ge",(spacing/2.0)+3*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge2","Ge",(spacing/2.0)+2*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge3","Ge",(spacing/2.0)+1*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge4","Ge",(spacing/2.0)+0*spacing);
          ConstructSNOLABDet(physicalWorld,"Si5","Si",-(spacing/2.0)-0*spacing);
          ConstructSNOLABDet(physicalWorld,"Si6","Si",-(spacing/2.0)-1*spacing);
          ConstructSNOLABDet(physicalWorld,"Si7","Si",-(spacing/2.0)-2*spacing);
          ConstructSNOLABDet(physicalWorld,"Si8","Si",-(spacing/2.0)-3*spacing);
        }
        else if(DesignNo==-11){
          //construct Diag9
	  //just a Ge/Si tower of SNOLAB dets 
	  //for cascade generation in the Si detectors ONLY
	  G4double spacing = 33.0+10.16; //use 0.4in spacing
	  //spacing negative switches order of towers
	  spacing *= -1;
          ConstructSNOLABDet(physicalWorld,"Ge1","Ge",(spacing/2.0)+3*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge2","Ge",(spacing/2.0)+2*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge3","Ge",(spacing/2.0)+1*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge4","Ge",(spacing/2.0)+0*spacing);
          ConstructSNOLABDet(physicalWorld,"Si5","Si",-(spacing/2.0)-0*spacing);
          ConstructSNOLABDet(physicalWorld,"Si6","Si",-(spacing/2.0)-1*spacing);
          ConstructSNOLABDet(physicalWorld,"Si7","Si",-(spacing/2.0)-2*spacing);
          ConstructSNOLABDet(physicalWorld,"Si8","Si",-(spacing/2.0)-3*spacing);
        }
        else if(DesignNo==-12){
          //construct Diag10
	  //just a Ge detector centered at zero
          ConstructGeDetSoudan(physicalWorld,0.0);
        }
        else if(DesignNo==-13){
	  //just a Si detector centered at zero
          ConstructR68SiDet(physicalWorld,0.0);
        }
        else if(DesignNo==-14){
          //construct Diag12
	  //just a Ne detector centered at zero
          ConstructNeDet(physicalWorld,0.0);
        }
        
	
	
	
	
	return physicalWorld;
}
void NeutReflect_DetectorConstruction::ConstructNULL(G4VPhysicalVolume *world)
{
        //keep running sum of downshift
        G4double downshift=0;

        DiskDisplacement=-6.35/2.0;
        //place the plastic
        ConstructSourcePlug(world,0.0,0.0,0.0);
        downshift+=6.35/2.0;

        //place the BeO disk
        downshift+=Bet/2.0;
        ConstructBePlug(world,0.0,0.0,-downshift);
        downshift+=Bet/2.0;

        return;
}
void NeutReflect_DetectorConstruction::ConstructDiag0(G4VPhysicalVolume *world)
{
        //make sure everything fits
        G4double llimit = BeR +10.0; //crude, it really depends on thickness too
        if(llimit>SphDetR)
          SphDetR = llimit;
        

        //add a spherical sensitive detector at origin
	G4Sphere* detSphere = new G4Sphere("detSphere_S",0.0,SphDetR,0.0,2*pi,0.0,pi);
	G4LogicalVolume* logicalDetSphere = new G4LogicalVolume(detSphere,Air,"detSphere_L",0,0,0);
	G4VPhysicalVolume* sphereDetWorld = new G4PVPlacement(negxrot,
								G4ThreeVector(0,0,0),
								"detSphere_P",
								logicalDetSphere,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttDetSphere = new G4VisAttributes(G4Colour(1.,0.,0.));
	VisAttDetSphere->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalDetSphere->SetVisAttributes(VisAttDetSphere);  
	// Make Invisible
	logicalDetSphere->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalDetSphere,"DetSphere");
        }

        //keep running sum of downshift
        G4double downshift=0;

        DiskDisplacement=-6.35/2.0;
        //place the plastic
        ConstructSourcePlug(sphereDetWorld,0.0,0.0,0.0);
        downshift+=6.35/2.0;

        //place the BeO disk
        downshift+=Bet/2.0;
        ConstructBePlug(sphereDetWorld,0.0,0.0,-downshift);
        downshift+=Bet/2.0;


        return;
}
void NeutReflect_DetectorConstruction::ConstructDiag1(G4VPhysicalVolume *world)
{
        //make sure everything fits
        G4double llimit = BeR +10.0; //crude, it really depends on thickness too
        if(llimit>SphDetR)
          SphDetR = llimit;
        

        //add a spherical sensitive detector at origin
	G4Sphere* detSphere = new G4Sphere("detSphere_S",0.0,SphDetR,0.0,2*pi,0.0,pi);
	G4LogicalVolume* logicalDetSphere = new G4LogicalVolume(detSphere,Air,"detSphere_L",0,0,0);
	G4VPhysicalVolume* sphereDetWorld = new G4PVPlacement(negxrot,
								G4ThreeVector(0,0,0),
								"detSphere_P",
								logicalDetSphere,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttDetSphere = new G4VisAttributes(G4Colour(1.,0.,0.));
	VisAttDetSphere->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalDetSphere->SetVisAttributes(VisAttDetSphere);  
	// Make Invisible
	logicalDetSphere->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalDetSphere,"DetSphere");
        }

        //keep running sum of downshift
        G4double downshift=0;

        DiskDisplacement=-5.2/2.0;
        //place the plastic
        ConstructStainlessPlug(sphereDetWorld,0.0,0.0,(17.6/2)-(5.2/2));
        downshift+=5.2/2.0;

        //place the BeO disk
        downshift+=Bet/2.0;
        ConstructBePlug(sphereDetWorld,0.0,0.0,-downshift);
        downshift+=Bet/2.0;


        return;
}
void NeutReflect_DetectorConstruction::ConstructDiag2(G4VPhysicalVolume *world)
{
        //make sure everything fits
        G4double llimit = BeR +10.0; //crude, it really depends on thickness too
        if(llimit>SphDetR)
          SphDetR = llimit;
       
	//make disk displacement
        DiskDisplacement=-6.35;

	//make a unit of inches
	static const G4double inch = 2.54*cm;

        //make a rotation to orient the system as other sources have
	G4ThreeVector r1 = G4ThreeVector(1.0,0,0);
	G4ThreeVector r2 = G4ThreeVector(0,-1,0);
	G4ThreeVector r3 = G4ThreeVector(0,0,-1);
	G4RotationMatrix *overallrot = new G4RotationMatrix(r1,r2,r3);
	overallrot->setRows(r1,r2,r3);

        G4ThreeVector disp = G4ThreeVector(0.0,-(-4.0+6*0.5+(0.25/2.0))*inch,-((0.25/2.0)+4.0)*inch).transform(*overallrot);
	
        //add a spherical sensitive detector at origin
	G4Sphere* detSphere = new G4Sphere("detSphere_S",0.0,SphDetR,0.0,2*pi,0.0,pi);
	G4LogicalVolume* logicalDetSphere = new G4LogicalVolume(detSphere,Air,"detSphere_L",0,0,0);
	G4VPhysicalVolume* sphereDetWorld = new G4PVPlacement(overallrot,
								disp + G4ThreeVector(0,0,0),
								"detSphere_P",
								logicalDetSphere,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttDetSphere = new G4VisAttributes(G4Colour(1.,0.,0.));
	VisAttDetSphere->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalDetSphere->SetVisAttributes(VisAttDetSphere);  
	// Make Invisible
	logicalDetSphere->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalDetSphere,"DetSphere");
        }

        //keep running sum of downshift
        G4double downshift=0;

        ConstructPhotoNSourceBox(sphereDetWorld,0.0,0.0,0.0);


        return;
}
void NeutReflect_DetectorConstruction::ConstructDiag3(G4VPhysicalVolume *world)
{
        //make sure everything fits
        G4double llimit = BeR +10.0; //crude, it really depends on thickness too
        if(llimit>SphDetR)
          SphDetR = llimit;
       
	//make disk displacement
        DiskDisplacement=-6.35;

	//make a unit of inches
	static const G4double inch = 2.54*cm;

        //make a rotation to orient the system as other sources have
	G4ThreeVector r1 = G4ThreeVector(1.0,0,0);
	G4ThreeVector r2 = G4ThreeVector(0,-1,0);
	G4ThreeVector r3 = G4ThreeVector(0,0,-1);
	G4RotationMatrix *overallrot = new G4RotationMatrix(r1,r2,r3);
	overallrot->setRows(r1,r2,r3);

        G4ThreeVector disp = G4ThreeVector(0.0,-(-4.0+6*0.5+(0.25/2.0))*inch,-((0.25/2.0)+4.0)*inch).transform(*overallrot);
	
        //add a spherical sensitive detector at origin
	G4Sphere* detSphere = new G4Sphere("detSphere_S",0.0,SphDetR,0.0,2*pi,0.0,pi);
	G4LogicalVolume* logicalDetSphere = new G4LogicalVolume(detSphere,Air,"detSphere_L",0,0,0);
	G4VPhysicalVolume* sphereDetWorld = new G4PVPlacement(overallrot,
								disp + G4ThreeVector(0,0,0),
								"detSphere_P",
								logicalDetSphere,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttDetSphere = new G4VisAttributes(G4Colour(1.,0.,0.));
	VisAttDetSphere->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalDetSphere->SetVisAttributes(VisAttDetSphere);  
	// Make Invisible
	logicalDetSphere->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalDetSphere,"DetSphere");
        }

        //keep running sum of downshift
        G4double downshift=0;

        ConstructPhotoNSourceBox(sphereDetWorld,0.0,0.0,0.0,"8020_medres.txt",0,1,false);


        return;
}
void NeutReflect_DetectorConstruction::ConstructDiag4(G4VPhysicalVolume *world)
{
        //make sure everything fits
        G4double llimit = BeR +10.0; //crude, it really depends on thickness too
        if(llimit>SphDetR)
          SphDetR = llimit;
       
	//make disk displacement
        DiskDisplacement=-6.35;

	//make a unit of inches
	static const G4double inch = 2.54*cm;

        //make a rotation to orient the system as other sources have
	G4ThreeVector r1 = G4ThreeVector(1.0,0,0);
	G4ThreeVector r2 = G4ThreeVector(0,-1,0);
	G4ThreeVector r3 = G4ThreeVector(0,0,-1);
	G4RotationMatrix *overallrot = new G4RotationMatrix(r1,r2,r3);
	overallrot->setRows(r1,r2,r3);

        G4ThreeVector disp = G4ThreeVector(0.0,-(-4.0+6*0.5+(0.25/2.0))*inch,-((0.25/2.0)+4.0)*inch).transform(*overallrot);
	
        //add a spherical sensitive detector at origin
	G4Sphere* detSphere = new G4Sphere("detSphere_S",0.0,SphDetR,0.0,2*pi,0.0,pi);
	G4LogicalVolume* logicalDetSphere = new G4LogicalVolume(detSphere,Air,"detSphere_L",0,0,0);
	G4VPhysicalVolume* sphereDetWorld = new G4PVPlacement(overallrot,
								disp + G4ThreeVector(0,0,0),
								"detSphere_P",
								logicalDetSphere,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttDetSphere = new G4VisAttributes(G4Colour(1.,0.,0.));
	VisAttDetSphere->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalDetSphere->SetVisAttributes(VisAttDetSphere);  
	// Make Invisible
	logicalDetSphere->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalDetSphere,"DetSphere");
        }

        //keep running sum of downshift
        G4double downshift=0;

        ConstructPhotoNSourceBox(sphereDetWorld,0.0,0.0,0.0,"8020_medres.txt",0,1,false,false);


        return;
}
void NeutReflect_DetectorConstruction::ConstructDesign0(G4VPhysicalVolume *world,G4double tl)
{
        //first, BeR can't be larger than inner vessel radius
        if(BeR>VesselR)
          BeR = VesselR;

        //keep running sum of downshift
        G4double downshift=0;

        DiskDisplacement=-6.35/2.0;
        //place the plastic
        ConstructSourcePlug(world,0.0,0.0,0.0);
        downshift+=6.35/2.0;

        //place the BeO disk
        downshift+=BeR/2.0;
        ConstructBePlug(world,0.0,0.0,-downshift);
        downshift+=BeR/2.0;

        //distance below origin to start
        G4double epsilon = downshift+1.0; //distance in mm

	//start by constructing lead slab
	G4Tubs* leadCylinder = new G4Tubs("PbCyl_S",0.0,76.2,tl*mfp_Pb_1_8/2,0.0,2*pi);

	G4LogicalVolume* logicalLeadCylinder = new G4LogicalVolume(leadCylinder,Pb,"PbCyl_L",0,0,0);

	G4VPhysicalVolume* cylinderLeadWorld = new G4PVPlacement(xrot,
								G4ThreeVector(0,-tl*mfp_Pb_1_8/2-epsilon,0),
								"PbCyl_P",
								logicalLeadCylinder,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttCyl = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	//VisAttCyl->SetForceSolid(false);
	VisAttCyl->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalLeadCylinder->SetVisAttributes(VisAttCyl);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        return;
}
void NeutReflect_DetectorConstruction::ConstructDesign1(G4VPhysicalVolume *world,G4double tl,
G4double tu,
G4double dr,
G4double g1,
G4double g2,
G4double hr,
G4double rr)
{
        //first, BeR can't be larger than inner vessel radius
        if(BeR>(76.2+dr+g2))
          BeR = 76.2+dr+g2;

        //keep running sum of downshift
        G4double downshift=0;

        DiskDisplacement=-6.35/2.0;
        //place the plastic
        ConstructSourcePlug(world,0.0,0.0,0.0);
        downshift+=6.35/2.0;

        //place the BeO disk
        downshift+=Bet/2.0;
        ConstructBePlug(world,0.0,0.0,-downshift);
        downshift+=Bet/2.0;

        //distance below origin to start
        G4double epsilon = downshift+1.0; //distance in mm

	//start by constructing lead slab
	G4Tubs* leadCylinder = new G4Tubs("PbCyl_S",0.0,76.2,tl*mfp_Pb_1_8/2,0.0,2*pi);

	G4LogicalVolume* logicalLeadCylinder = new G4LogicalVolume(leadCylinder,Pb,"PbCyl_L",0,0,0);

	G4VPhysicalVolume* cylinderLeadWorld = new G4PVPlacement(xrot,
								G4ThreeVector(0,-tl*mfp_Pb_1_8/2-epsilon,0),
								"PbCyl_P",
								logicalLeadCylinder,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttCyl = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	//VisAttCyl->SetForceSolid(false);
	VisAttCyl->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalLeadCylinder->SetVisAttributes(VisAttCyl);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        //add another slab above
	G4Tubs* leadCylinder2 = new G4Tubs("PbCyl2_S",0.0,76.2+dr,tu/2,0.0,2*pi);

	G4LogicalVolume* logicalLeadCylinder2 = new G4LogicalVolume(leadCylinder2,Pb,"PbCyl_L",0,0,0);

	G4VPhysicalVolume* cylinder2LeadWorld = new G4PVPlacement(xrot, 
								G4ThreeVector(0,+tu/2+g1,0),
								"PbCyl2_P",
								logicalLeadCylinder2,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttCyl2 = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	//VisAttCyl2->SetForceSolid(false);
	VisAttCyl2->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalLeadCylinder2->SetVisAttributes(VisAttCyl2);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);
	
	//enclosing slab 
	G4Tubs* leadCylinder3 = new G4Tubs("PbCyl3_S",76.2+dr+g2,76.2+dr+g2+rr,hr/2,0.0,2*pi);

	G4LogicalVolume* logicalLeadCylinder3 = new G4LogicalVolume(leadCylinder3,Pb,"PbCyl_L",0,0,0);

	G4VPhysicalVolume* cylinder3LeadWorld = new G4PVPlacement(xrot, 
								G4ThreeVector(0,0,0),
								"PbCyl3_P",
								logicalLeadCylinder3,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttCyl3 = new G4VisAttributes(G4Colour(0.,0.,1.));
	//VisAttCyl2->SetForceSolid(false);
	VisAttCyl3->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalLeadCylinder3->SetVisAttributes(VisAttCyl3);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible); 
	
        return;
}
void NeutReflect_DetectorConstruction::ConstructDesign2(G4VPhysicalVolume *world,G4double tl,
G4double gedist,
G4double gerad,
G4double gethk,
G4double beorad)
{
        //first, BeR can't be larger than inner vessel radius
        if(BeR>VesselR)
          BeR = VesselR;

        //keep running sum of downshift
        G4double downshift=0;

        DiskDisplacement=-6.35/2.0;
        //place the plastic
        ConstructSourcePlug(world,0.0,0.0,0.0);
        downshift+=6.35/2.0;

        //place the BeO disk
        downshift+=Bet/2.0;
        ConstructBePlug(world,0.0,0.0,-downshift);
        downshift+=Bet/2.0;

        //distance below origin to start
        G4double epsilon = downshift+1.0; //distance in mm

	//start by constructing lead slab cone
        G4double leadlargerad = (tl*mfp_Pb_1_8)*gerad*(1/(gedist-(gethk/2)));
 	G4cout << "Lead large rad is: " << leadlargerad << G4endl;
	G4Cons* leadCone = new G4Cons("PbCone_S",0.0,leadlargerad,0.0,5.0,tl*mfp_Pb_1_8/2,0.0,2*pi);

	G4LogicalVolume* logicalLeadCone = new G4LogicalVolume(leadCone,Pb,"PbCone_L",0,0,0);

	G4VPhysicalVolume* coneLeadWorld = new G4PVPlacement(xrot,
								G4ThreeVector(0,-tl*mfp_Pb_1_8/2-epsilon,0),
								"PbCone_P",
								logicalLeadCone,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttCone = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	//VisAttCone->SetForceSolid(false);
	VisAttCone->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalLeadCone->SetVisAttributes(VisAttCone);  
	// Make Invisible
	//logicalCone->SetVisAttributes(G4VisAttributes::Invisible);

        return;
}
void NeutReflect_DetectorConstruction::ConstructDesign3(G4VPhysicalVolume *world,G4double tl,
G4double dpipe,
G4double tsteelplate,
G4double tsrcholder, 
G4double lBeO,
G4double ldet,
G4double ddet)
{
        //keep some other hard-coded params
        G4double botplatemargin=30.0,topplatemargin=60.0,pipethk=3.175,drod=25.0,botplaterodmargin=5.0; 
        G4double tsrc=6.35; 

        //keep running sum of downshift
        G4double downshift=0;

        DiskDisplacement=-6.35/2.0;
        //place the plastic
        ConstructSourcePlug(world,0.0,0.0,0.0);
        downshift+=6.35/2.0;

        //place the BeO disk
        downshift+=Bet/2.0;
        //ConstructBePlug(world,0.0,0.0,-downshift);
        downshift+=Bet/2.0;

        //distance below origin to start
        G4double epsilon = downshift+1.0; //distance in mm

        // define some useful geometrical objects
	G4ThreeVector shift = G4ThreeVector(0,0,0); //for making sphere-cylinder combination
	G4RotationMatrix* rm = new G4RotationMatrix(); // rotation of second object in combination
        rm->setRows(G4ThreeVector(1,0,0),G4ThreeVector(0,1,0),G4ThreeVector(0,0,1));

        //--------------------------------------------------
        // create the framework for the source holder and place it
        //--------------------------------------------------
	
        //start by constructing steel bottom plate
        G4double shiftforbottomplate = (tl*mfp_Pb_1_8) + tsrcholder;
	G4Tubs* steelCylinderBot = new G4Tubs("SteelCylBot_S",dpipe/2.0,tl*mfp_Pb_1_8 + botplatemargin ,(tsteelplate)/2.0,0.0,2*pi);
	//G4Tubs* steelCylinderBot = new G4Tubs("SteelCylBot_S",0.0,tl*mfp_Pb_1_8 + 30.0 ,(tsteelplate)/2.0,0.0,2*pi);

        //make a pipe
	G4Tubs* steelCylinderPipe = new G4Tubs("SteelCylPipe_S",dpipe/2.0,(dpipe/2.0)+pipethk,(tl*mfp_Pb_1_8)/2.0,0.0,2*pi);
	//G4Tubs* steelCylinderPipe = new G4Tubs("SteelCylPipe_S",0.0,(dpipe/2.0)+3.175,(tl*mfp_Pb_1_8)/2.0,0.0,2*pi);

        //make a rod
	G4Tubs* steelCylinderRod = new G4Tubs("SteelCylRod_S",0.0,drod/2.0,(2*tl*mfp_Pb_1_8 +2*tsrcholder)/2.0,0.0,2*pi);

        //make a top plate
	G4Tubs* steelCylinderTop = new G4Tubs("SteelCylTop_S",0.0,tl*mfp_Pb_1_8 + topplatemargin,tsteelplate/2.0,0.0,2*pi);

        //join bottom plate and pipe
        shift = G4ThreeVector(0.0,0.0,-((tl*mfp_Pb_1_8)/2.0)-(tsteelplate/2.0));
	G4Transform3D transform(*rm,shift);
        G4VSolid *botPlatePipe = new G4UnionSolid("botPlatePipe",steelCylinderBot,steelCylinderPipe,transform);
        //G4cout << "BLAHHHH DISTANCE: " << botPlatePipe->DistanceToIn(G4ThreeVector(0,0,300.0))<<":"<<steelCylinderBot->DistanceToIn(G4ThreeVector(0,0,300.0)) << G4endl;
        //G4cout << "BLAHHHH DISTANCE: " << ((tl*mfp_Pb_1_8)+tsteelplate)/2.0 << G4endl;

        //join the supports in a symmetric pattern
        G4int nsupp=6; //make it even
        G4double angle=0;
        G4double radius=(tl*mfp_Pb_1_8)+botplatemargin - (drod/2.0) - botplaterodmargin;
        G4VSolid *botPlatePipeandRods = botPlatePipe;
        for(G4int i=0;i<nsupp;i++){
          shift = G4ThreeVector(radius*cos(angle),radius*sin(angle),-(2*tl*mfp_Pb_1_8 + 2*tsrcholder)/2.0 - (tsteelplate/2.0));
	  G4Transform3D transform(*rm,shift);
          G4VSolid *savePtr = botPlatePipeandRods;
          botPlatePipeandRods = new G4UnionSolid("botPlatePipeandRods",botPlatePipeandRods,steelCylinderRod,transform);
          angle += (2*pi)/nsupp;

        }

        //join the top plate
        shift = G4ThreeVector(0.0,0.0,-((2*tl*mfp_Pb_1_8 + 2*tsrcholder))-(tsteelplate/2.0)-(tsteelplate/2.0));
	G4Transform3D transform1(*rm,shift);
        G4VSolid *botPlatePipeandRodsandTop = new G4UnionSolid("botPlatePipeandRodsandTop",botPlatePipeandRods,steelCylinderTop,transform1);


        //place the assembly 
	G4LogicalVolume* logicalSteelFrame = new G4LogicalVolume(botPlatePipeandRodsandTop,Steel,"SteelFrame_L",0,0,0);

	G4VPhysicalVolume* frameSteelWorld = new G4PVPlacement(xrot,
								G4ThreeVector(0,-(tsteelplate/2.0)-(tl*mfp_Pb_1_8)-(tsrcholder),0),
								"SteelFrame_P",
								logicalSteelFrame,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttFrame = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	//VisAttCone->SetForceSolid(false);
	VisAttFrame->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalSteelFrame->SetVisAttributes(VisAttFrame);  
	// Make Invisible
	//logicalCone->SetVisAttributes(G4VisAttributes::Invisible);


        //--------------------------------------------------
        // create the source holder and place it
        //--------------------------------------------------

        //start by constructing steel bottom plate
	G4Tubs* steelSrcPlatterBot = new G4Tubs("SteelSrcPlatterBot_S",dpipe/2.0,tl*mfp_Pb_1_8 ,(tsrcholder/2.0)/2.0,0.0,2*pi);

        //constructing steel top plate
	G4Tubs* steelSrcPlatterTop = new G4Tubs("SteelSrcPlatterTop_S",0.0,tl*mfp_Pb_1_8 ,(tsrcholder/2.0)/2.0,0.0,2*pi);

        //constructing lead innards 
	G4Tubs* leadSrcPlatterIn = new G4Tubs("LeadSrcPlatterIn_S",12.7,tl*mfp_Pb_1_8 ,(tsrcholder)/2.0,0.0,2*pi);

        //place the source holder 
	G4LogicalVolume* logicalSrcPlatterBot = new G4LogicalVolume(steelSrcPlatterBot,Steel,"SteelSrcPlatterBot_L",0,0,0);

	G4VPhysicalVolume* srcPlatterBotSteelWorld = new G4PVPlacement(xrot,
								G4ThreeVector(0,-tsrcholder*(3.0/4.0),0),
								"SteelSrcPlatterBot_P",
								logicalSrcPlatterBot,
								world,
								false,
								0);

	logicalSrcPlatterBot->SetVisAttributes(VisAttFrame);  

	G4LogicalVolume* logicalSrcPlatterTop = new G4LogicalVolume(steelSrcPlatterTop,Steel,"SteelSrcPlatterTop_L",0,0,0);

	G4VPhysicalVolume* srcPlatterTopSteelWorld = new G4PVPlacement(xrot,
								G4ThreeVector(0,+tsrcholder*(3.0/4.0),0),
								"SteelSrcPlatterTop_P",
								logicalSrcPlatterTop,
								world,
								false,
								0);

	logicalSrcPlatterTop->SetVisAttributes(VisAttFrame);  

	G4LogicalVolume* logicalSrcPlatterIn = new G4LogicalVolume(leadSrcPlatterIn,Pb,"LeadSrcPlatterTop_L",0,0,0);

	G4VPhysicalVolume* srcPlatterInPbWorld = new G4PVPlacement(xrot,
								G4ThreeVector(0,0,0),
								"LeadSrcPlatterIn_P",
								logicalSrcPlatterIn,
								world,
								false,
								0);

	//logicalSrcPlatterIn->SetVisAttributes(VisAttFrame);  
	// Visualization attributes
	G4VisAttributes* VisAttIn = new G4VisAttributes(G4Colour(0.7,0.0,0.3));
	//VisAttCone->SetForceSolid(false);
	VisAttIn->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalSrcPlatterIn->SetVisAttributes(VisAttIn);  
	// Make Invisible
	//logicalCone->SetVisAttributes(G4VisAttributes::Invisible);

        //--------------------------------------------------
        // create the shielding and place it
        //--------------------------------------------------

        //constructing lead bottom 
	G4Tubs* leadShieldBot = new G4Tubs("LeadShieldBot_S",(dpipe/2.0)+pipethk,tl*mfp_Pb_1_8 ,(tl*mfp_Pb_1_8)/2.0,0.0,2*pi);

        //constructing lead top 
	G4Tubs* leadShieldTop = new G4Tubs("LeadShieldTop_S",0.0,tl*mfp_Pb_1_8 ,(tl*mfp_Pb_1_8)/2.0,0.0,2*pi);

        //place the shields 
	G4LogicalVolume* logicalShieldBot = new G4LogicalVolume(leadShieldBot,Pb,"LeadShieldBot_L",0,0,0);

	G4VPhysicalVolume* shieldBotLeadWorld = new G4PVPlacement(xrot,
								G4ThreeVector(0,-tsrcholder-(tl*mfp_Pb_1_8/2.0),0),
								"LeadShieldBot_P",
								logicalShieldBot,
								world,
								false,
								0);

	logicalShieldBot->SetVisAttributes(VisAttFrame);  

	G4LogicalVolume* logicalShieldTop = new G4LogicalVolume(leadShieldTop,Pb,"LeadShieldTop_L",0,0,0);

	G4VPhysicalVolume* shieldTopLeadWorld = new G4PVPlacement(xrot,
								G4ThreeVector(0,+tsrcholder+(tl*mfp_Pb_1_8/2.0),0),
								"LeadShieldTop_P",
								logicalShieldTop,
								world,
								false,
								0);

	logicalShieldTop->SetVisAttributes(VisAttFrame);  

        return;
}
void NeutReflect_DetectorConstruction::ConstructSourcePlug(G4VPhysicalVolume *world,G4double xpos,
G4double ypos,
G4double zpos)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        G4cout << world->GetName() << G4endl;
        if(world->GetName()=="world_P"){
          thisrot=xrot;
        }
        else{
          thisrot=nullrot;
        }

        //displacement vector
        G4ThreeVector disp = G4ThreeVector(xpos,ypos,zpos).transform(*thisrot);

        //create a disk
	G4Tubs* srcCylinder = new G4Tubs("srcCyl_S",0.0,12.7,3.175,0.0,2*pi);
	G4LogicalVolume* logicalSrcCylinder = new G4LogicalVolume(srcCylinder,Plastic,"SrcCyl_L",0,0,0);
	G4VPhysicalVolume* cylinderSrcWorld = new G4PVPlacement(thisrot, 
								disp,
								"SrcCyl_P",
								logicalSrcCylinder,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttSrcCyl = new G4VisAttributes(G4Colour(7.,0.,3.));
	VisAttSrcCyl->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalSrcCylinder->SetVisAttributes(VisAttSrcCyl);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        return;
}
void NeutReflect_DetectorConstruction::ConstructStainlessPlug(G4VPhysicalVolume *world,G4double xpos,
G4double ypos,
G4double zpos)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        G4cout << world->GetName() << G4endl;
        if(world->GetName()=="world_P"){
          thisrot=xrot;
        }
        else{
          thisrot=nullrot;
        }

        //displacement vector
        G4ThreeVector disp = G4ThreeVector(xpos,ypos,zpos).transform(*thisrot);

        //create a disk
	G4Tubs* srcCylinder = new G4Tubs("srcCyl_S",0.0,3.2,17.6/2,0.0,2*pi);
	G4LogicalVolume* logicalSrcCylinder = new G4LogicalVolume(srcCylinder,stainlessSteel,"SrcCyl_L",0,0,0);
	G4VPhysicalVolume* cylinderSrcWorld = new G4PVPlacement(thisrot, 
								disp,
								"SrcCyl_P",
								logicalSrcCylinder,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttSrcCyl = new G4VisAttributes(G4Colour(140/255.0,253/255.0,153/255.0));
	VisAttSrcCyl->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalSrcCylinder->SetVisAttributes(VisAttSrcCyl);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        return;
}
void NeutReflect_DetectorConstruction::ConstructBePlug(G4VPhysicalVolume *world,G4double xpos,
G4double ypos,
G4double zpos)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        if(world->GetName()=="world_P")
          thisrot=xrot;
        else
          thisrot=nullrot;

        //displacement vector
        G4ThreeVector disp = G4ThreeVector(xpos,ypos,zpos).transform(*thisrot);
        //G4ThreeVector disp = G4ThreeVector(xpos,ypos,zpos);
        //G4cout << disp << G4endl;
        disp = G4ThreeVector(xpos,ypos,zpos);
        //G4cout << disp << G4endl;

        //create a disk (software-coded radius and thickness BeR and Bet)
	G4Tubs* beCylinder = new G4Tubs("beCyl_S",0.0,BeR,Bet/2.0,0.0,2*pi);
	G4LogicalVolume* logicalBeCylinder;
        if(!UseBePure)
          logicalBeCylinder = new G4LogicalVolume(beCylinder,BeO,"BeCyl_L",0,0,0);
        else
          logicalBeCylinder = new G4LogicalVolume(beCylinder,BePure,"BeCyl_L",0,0,0);
	G4VPhysicalVolume* cylinderBeWorld = new G4PVPlacement(thisrot, 
								disp,
								"BeCyl_P",
								logicalBeCylinder,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttBeCyl = new G4VisAttributes(G4Colour(7.,3.,0.));
	VisAttBeCyl->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalBeCylinder->SetVisAttributes(VisAttBeCyl);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        return;
}
void NeutReflect_DetectorConstruction::ConstructGeDet(G4VPhysicalVolume *world,G4double dist)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        if(world->GetName()=="world_P")
          thisrot=xrot;
        else
          thisrot=nullrot;

	//make a germanium disk R=3" thickness=1.5"
	//FIXME hardcoded size
	G4Tubs* geCylinder = new G4Tubs("geCyl_S",0.0,76.2,(25.4+12.7)/2,0.0,2*pi);
	G4LogicalVolume* logicalGeCylinder = new G4LogicalVolume(geCylinder,Ge,"GeCyl_L",0,0,0);
	G4VPhysicalVolume* cylinderGeWorld = new G4PVPlacement(thisrot,
								G4ThreeVector(0,-dist,0),
								"GeCyl_P",
								logicalGeCylinder,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttGeCyl = new G4VisAttributes(G4Colour(1.,0.,0.));
	//VisAttCyl2->SetForceSolid(false);
	VisAttGeCyl->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalGeCylinder->SetVisAttributes(VisAttGeCyl);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalGeCylinder,"GeDet");
        }
        return;
}
void NeutReflect_DetectorConstruction::ConstructGeDetSoudan(G4VPhysicalVolume *world,G4double dist)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        if(world->GetName()=="world_P")
          thisrot=xrot;
        else
          thisrot=nullrot;

	//make a germanium disk R=1.5" thickness=1"
	//FIXME hardcoded size
	G4Tubs* geCylinder = new G4Tubs("geCyl_S",0.0,76.2/2.0,(25.4)/2,0.0,2*pi);
	G4LogicalVolume* logicalGeCylinder = new G4LogicalVolume(geCylinder,Ge,"GeCyl_L",0,0,0);
	G4VPhysicalVolume* cylinderGeWorld = new G4PVPlacement(thisrot,
								G4ThreeVector(0,-dist,0),
								"GeCyl_P",
								logicalGeCylinder,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttGeCyl = new G4VisAttributes(G4Colour(1.,0.,0.));
	//VisAttCyl2->SetForceSolid(false);
	VisAttGeCyl->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalGeCylinder->SetVisAttributes(VisAttGeCyl);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalGeCylinder,"GeDet");
        }
        return;
}
void NeutReflect_DetectorConstruction::ConstructSiDet(G4VPhysicalVolume *world,G4double dist)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        if(world->GetName()=="world_P")
          thisrot=xrot;
        else
          thisrot=nullrot;

	//make a silicon disk R=3" thickness=1.5"
	//FIXME hardcoded size
	G4Tubs* siCylinder = new G4Tubs("siCyl_S",0.0,76.2,(25.4+12.7)/2,0.0,2*pi);
	G4LogicalVolume* logicalSiCylinder = new G4LogicalVolume(siCylinder,Si,"SiCyl_L",0,0,0);
	G4VPhysicalVolume* cylinderGeWorld = new G4PVPlacement(thisrot,
								G4ThreeVector(0,-dist,0),
								"SiCyl_P",
								logicalSiCylinder,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttSiCyl = new G4VisAttributes(G4Colour(1.,1.,0.));
	//VisAttCyl2->SetForceSolid(false);
	VisAttSiCyl->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalSiCylinder->SetVisAttributes(VisAttSiCyl);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalSiCylinder,"SiDet");
        }
        return;
}
void NeutReflect_DetectorConstruction::ConstructR68SiDet(G4VPhysicalVolume *world,G4double dist)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        if(world->GetName()=="world_P")
          thisrot=xrot;
        else
          thisrot=nullrot;

	//make a silicon disk R=3" thickness=1.5"
	//FIXME hardcoded size
	G4Tubs* siCylinder = new G4Tubs("siCyl_S",0.0,50.0,(33.0)/2,0.0,2*pi);
	G4LogicalVolume* logicalSiCylinder = new G4LogicalVolume(siCylinder,Si,"SiCyl_L",0,0,0);
	G4VPhysicalVolume* cylinderGeWorld = new G4PVPlacement(thisrot,
								G4ThreeVector(0,-dist,0),
								"SiCyl_P",
								logicalSiCylinder,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttSiCyl = new G4VisAttributes(G4Colour(1.,1.,0.));
	//VisAttCyl2->SetForceSolid(false);
	VisAttSiCyl->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalSiCylinder->SetVisAttributes(VisAttSiCyl);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalSiCylinder,"R68SiDet");
        }
        return;
}
void NeutReflect_DetectorConstruction::ConstructSNOLABDet(G4VPhysicalVolume *world,G4String name, G4String mat,G4double dist)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        if(world->GetName()=="world_P")
          thisrot=xrot;
        else
          thisrot=nullrot;

	//make some names for volumes
	G4String shapename = name+"_S";
	G4String lvolname = name+"_L";
	G4String pvolname = name+"_P";

	//make a germanium disk R=3" thickness=1.5"
	//FIXME hardcoded size
	G4Tubs* Cylinder = new G4Tubs(shapename,0.0,(100.0/2.0),(33.0)/2,0.0,2*pi);
	G4LogicalVolume* logicalCylinder;
        if(mat=="Ge")	
	  logicalCylinder = new G4LogicalVolume(Cylinder,Ge,lvolname,0,0,0);
	else if(mat=="Si")	
	  logicalCylinder = new G4LogicalVolume(Cylinder,Si,lvolname,0,0,0);
	else //use Ge as default?	
	  logicalCylinder = new G4LogicalVolume(Cylinder,Ge,lvolname,0,0,0);

	G4VPhysicalVolume* cylinderWorld = new G4PVPlacement(thisrot,
								G4ThreeVector(0,-dist,0),
								pvolname,
								logicalCylinder,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttCyl;
	if(mat=="Ge")
          VisAttCyl = new G4VisAttributes(G4Colour(1.,0.,0.));
	else if(mat=="Si")
          VisAttCyl = new G4VisAttributes(G4Colour(1.,1.,0.));
	else
          VisAttCyl = new G4VisAttributes(G4Colour(1.,0.,0.));
	VisAttCyl->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalCylinder->SetVisAttributes(VisAttCyl);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalCylinder,name);
        }
        return;
}
void NeutReflect_DetectorConstruction::ConstructArDet(G4VPhysicalVolume *world,G4double dist)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        if(world->GetName()=="world_P")
          thisrot=xrot;
        else
          thisrot=nullrot;

	//make an argon sphere of 1 m at STP
	//FIXME hardcoded size
	G4Sphere* arSphere = new G4Sphere("arSphere_S",0.0,1000.0,0.0,2*pi,0.0,pi);
	G4LogicalVolume* logicalArSphere = new G4LogicalVolume(arSphere,Argon,"ArSph_L",0,0,0);
	G4VPhysicalVolume* cylinderArWorld = new G4PVPlacement(thisrot,
								G4ThreeVector(0,-dist,0),
								"ArSph_P",
								logicalArSphere,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttArSph = new G4VisAttributes(G4Colour(255./255.,127./255.,80./255.));
	VisAttArSph->SetForceSolid(false);
	//VisAttArSph->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalArSphere->SetVisAttributes(VisAttArSph);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalArSphere,"ArDet");
        }
        return;
}
void NeutReflect_DetectorConstruction::ConstructNeDet(G4VPhysicalVolume *world,G4double dist)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        if(world->GetName()=="world_P")
          thisrot=xrot;
        else
          thisrot=nullrot;

	//make an argon sphere of 1 m at STP
	//FIXME hardcoded size
	G4Sphere* neSphere = new G4Sphere("neSphere_S",0.0,1000.0,0.0,2*pi,0.0,pi);
	G4LogicalVolume* logicalNeSphere = new G4LogicalVolume(neSphere,Neon,"NeSph_L",0,0,0);
	G4VPhysicalVolume* cylinderNeWorld = new G4PVPlacement(thisrot,
								G4ThreeVector(0,-dist,0),
								"NeSph_P",
								logicalNeSphere,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttNeSph = new G4VisAttributes(G4Colour(255./255.,127./255.,80./255.));
	VisAttNeSph->SetForceSolid(false);
	//VisAttNeSph->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalNeSphere->SetVisAttributes(VisAttNeSph);  
	// Make Invisible
	//logicalCylinder->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>0){
          ConstructGenericSensitive(logicalNeSphere,"NeDet");
        }
        return;
}
void NeutReflect_DetectorConstruction::ConstructSimpleVessel(G4VPhysicalVolume *world,G4double vthk,
G4double vrad,
G4double vh)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        if(world->GetName()=="world_P")
          thisrot=xrot;
        else
          thisrot=nullrot;

	//make a steel
	G4Tubs* steelCylinderOut = new G4Tubs("stCylOut_S",0.0,vrad+vthk,vh/2,0.0,2*pi);
	G4Tubs* steelCylinderIn = new G4Tubs("stCylIn_S",0.0,vrad,(vh-vthk*2)/2,0.0,2*pi);
	G4VSolid* steelShell = new G4SubtractionSolid("stShell_S",steelCylinderOut,
                                                                        steelCylinderIn,
                                                                        0,
                                                                        G4ThreeVector(0,0,0));
	
	G4LogicalVolume* logicalSteelShell = new G4LogicalVolume(steelShell,stainlessSteel,"stShell_L",0,0,0);
	G4VPhysicalVolume* shellSteelWorld = new G4PVPlacement(thisrot, //no rotation
			     				       G4ThreeVector(0,0,0),
							       "stShell_P",
								logicalSteelShell,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttStShell = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
	//VisAttStShell->SetForceSolid(false);
	VisAttStShell->SetForceWireframe(true);  //I want a Wireframe of the me
	logicalSteelShell->SetVisAttributes(VisAttStShell);  
	// Make Invisible
	//logicalSteelShell->SetVisAttributes(G4VisAttributes::Invisible);

        //make this sensitive
        if(ConstructGenericSensitiveInt>1){
          ConstructGenericSensitive(logicalSteelShell,"SteelShell");
        }
        return;

}
void NeutReflect_DetectorConstruction::ConstructPhotoNSourceBox(G4VPhysicalVolume *world,G4double xpos,
G4double ypos,
G4double zpos,
G4String file8020,
G4int bshelf,
G4int gshelf,
G4bool dolead,
G4bool reston)
{
        /* ****************************************************************************************************
	 * There are many components here, they should probably be listed and taken care to implement each one
	 * correctly.  The level of detail here might be overkill, but then we are interested in knowing very
	 * precisely how low energy (100's of keV) neutrons scatter in and around the source.  Constructing the
	 * piece is based off of the "Bottom Plate," everything is placed relative to that.
	 *
	 *   Bottom Plate (aluminum)
	 *   Side Plate (aluminum) (2X)
	 *   Bridge Coupling Plate (aluminum)
	 *   Bridge Structure (aluminum)
	 *   Gamma Strut (aluminum)
	 *   Gamma Bar  (aluminum)
	 *   Be Strut  (aluminum)
	 *   Be Source Plate (aluminum)
	 *   Be Source Wafer (beryllium)
	 *   Source Capsule (stainless steel)
	 *   
	 *   NOTE: One can turn off the piece or just the visualization of that piece, these are NOT the same
	 *   as turning off the piece will have it removed from the simulation entirely but turning off the 
	 *   visualization will simply make it invisible to all geometry viewers (i.e. particles are still 
	 *   tracked through it). 
	 *
	 *   NOTE: (body coordinate system)
	 *
	 *   NOTE: typical units are inches for the design, I will define this unit and use it explicitly
	 *   throughout the construction.
	 *
	 */

	//make a unit of inches
	static const G4double inch = 2.54*cm;

	//turn on and off certain pieces
	G4bool bpOn=reston; //Bottom Plate
        G4bool spOn=reston; //Side Plate(s)
	G4bool bcpOn=reston; //Bridge Coupling Plate
	G4bool bstructOn=reston; //Bridge Structure 
	G4bool gsOn=reston; //Gamma Strut
	G4bool gbOn=reston; //Gamma Bar
	G4bool besOn=reston; //Be Strut
	G4bool bespOn=reston; //Be Source Plate
	G4bool bswOn=true; //Be Source Wafer
	G4bool scOn=reston; //Source Capsule
	G4bool lbrickOn=dolead; //Lead Bricks

	//turn on and off certain visualizations 
	G4bool bpVis=false; //Bottom Plate
        G4bool spVis=true; //Side Plate(s)
	G4bool bcpVis=true; //Bridge Coupling Plate
	G4bool bstructVis=true; //Bridge Structure 
	G4bool gsVis=true; //Gamma Strut
	G4bool gbVis=true; //Gamma Bar
	G4bool besVis=true; //Be Strut
	G4bool bespVis=true; //Be Source Plate
	G4bool bswVis=true; //Be Source Wafer
	G4bool scVis=true; //Source Capsule
	G4bool lbrickVis=true; //Lead Bricks

	//make a switch for doing things more fancy
	G4bool doFancy=true;

        //make a rotation to orient the system as other sources have
	G4ThreeVector r1 = G4ThreeVector(1.0,0,0);
	G4ThreeVector r2 = G4ThreeVector(0,0.0,1.0);
	G4ThreeVector r3 = G4ThreeVector(0,-1.0,0.0);
	G4RotationMatrix *overallrot = new G4RotationMatrix(r1,r2,r3);
	overallrot->setRows(r1,r2,r3);

        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        G4cout << world->GetName() << G4endl;
        if(world->GetName()=="world_P"){
          thisrot=xrot;
        }
        else{
          thisrot=nullrot;
        }

	xrot->print(G4cout);
	thisrot->print(G4cout);
        //displacement vector
        G4ThreeVector disp = G4ThreeVector(xpos,ypos,zpos).transform(*thisrot);
        //G4ThreeVector disp = G4ThreeVector(0.0,-(-4.0+6*0.5+(0.25/2.0))*inch,-((0.25/2.0)+4.0)*inch).transform(*overallrot);
        //G4ThreeVector disp_unrot = G4ThreeVector(0.0,-(-4.0+6*0.5+(0.25/2.0))*inch,-((0.25/2.0)+4.0)*inch);

        //various geometry objects needed for computations
        G4ThreeVector shift,row1,row2,row3;

        /****************************Set up Bottom Plate Construction*******************************************/
        //bottom plate subtractions
        G4Tubs* small_hole = new G4Tubs("small_hole",0.0,4.32/2.0,10,0.0,2*pi);
        G4Tubs* large_hole = new G4Tubs("large_hole",0.0,4.32/2.0,10,0.0,2*pi);
        G4Box* shelf = new G4Box("shelf",209.55/2.0,6.35/2.0,10.0/2.0);

        //bottom plate
	G4VSolid* botPlate = new G4Box("botPlate_S",228.6/2.0,203.2/2.0,6.35/2.0);

        //make holes around
        float corner_x=-114.3,corner_y=-101.6;
        float ystep=25.4;
        float yridgestep=12.7;
        float x1=4.76,x2=223.8;
        shift = G4ThreeVector(0.0,0.0,0.0);
        for(int i=0;i<7;i++){

           //get the iteration
      	   std::ostringstream oss;
	   oss << i;
        
           //small holes on the "left"
           shift.setX(corner_x+x1);
           shift.setY(corner_y+(i+1)*ystep);
           botPlate = new G4SubtractionSolid("botPlate_S_left_"+oss.str(),botPlate,small_hole,NULL,shift);

           //small holes on the "right"
           shift.setX(corner_x+x2);
           shift.setY(corner_y+(i+1)*ystep);
           botPlate = new G4SubtractionSolid("botPlate_S_right_"+oss.str(),botPlate,small_hole,NULL,shift);
        }
        for(int i=0;i<8;i++){
        
           //get the iteration
	   std::ostringstream oss;
	   oss << i;

           //ridges
           shift.setX(0.0);
           shift.setY(corner_y+(i+1)*yridgestep);
           shift.setZ((10.0/2.0)+(6.35/2.0)-1.59); //1/16" depth
           //shift.setZ((10.0/2.0)-1.59); // deeper for visualization ease
           botPlate = (G4VSolid*) new G4SubtractionSolid("botPlate_S_ridge_"+oss.str(),botPlate,shelf,NULL,shift);

        }

        /****************************End Bottom Plate Construction*******************************************/

        /****************************Set up Side Plate Construction******************************************/
        //side plate subtractions
        G4Tubs* small_side_hole = new G4Tubs("small_side_hole",0.0,3.45/2.0,19.05,0.0,2*pi);

        //side plate
	G4VSolid* sidePlate = new G4Box("sidePlate_S",203.2/2.0,203.2/2.0,9.525/2.0);
      
        //need a rotation for these holes
        G4RotationMatrix *holerot;
	row1 = G4ThreeVector(0,0,1);
	row2 = G4ThreeVector(0,1,0);
	row3 = G4ThreeVector(-1,0,0);
	holerot = new G4RotationMatrix(row1,row2,row3);
	holerot->setRows(row1,row2,row3);

        //make holes around
        float s_corner_x=-101.6,s_corner_y=-101.6;
        float s_ystep=25.4;
        shift = G4ThreeVector(0.0,0.0,0.0);
        for(int i=0;i<7;i++){
        
           //get the iteration
	   std::ostringstream oss;
	   oss << i;

           //small holes on the "left"
           shift.setX(s_corner_x);
           shift.setY(s_corner_y+(i+1)*s_ystep);
           sidePlate = new G4SubtractionSolid("sidePlate_S_left_"+oss.str(),sidePlate,small_side_hole,holerot,shift);

        }
        /****************************End Side Plate Construction*******************************************/

        /****************************Set up Bridge Coupling Plate Construction*****************************/
        //bc plate subtractions
        G4Box* large_bcp_hole = new G4Box("large_bcp_hole",(3.0*inch/2.0),(3.0*inch/2.0),25.4);

        //bc plate
	G4VSolid* bcPlate = new G4Box("bcPlate_S",292.1/2.0,241.3/2.0,6.35/2.0);
      
        //make holes around
        //float s_ystep=25.4;
        shift = G4ThreeVector(0.0,0.0,0.0);
        bcPlate = new G4SubtractionSolid("bcPlate_S",bcPlate,large_bcp_hole,NULL,shift);
        /****************************End Bridge Coupling Plate Construction*********************************/

        /****************************Set up Gamma Strut Construction****************************************/
        //gamma strut subtractions
        G4Tubs* gammastrut_through_hole = new G4Tubs("gammastrut_through_hole",0.0,4.318/2.0,25.4,0.0,2*pi);
        G4Tubs* gammastrut_center_hole = new G4Tubs("gammastrut_center_hole",0.0,6.527/2.0,25.4,0.0,2*pi);
	G4Box* gammastrut_cutout = new G4Box("gammastrut_cutout",12.7/2.0,25.4,12.7/2.0);

        //gamma strut
	G4VSolid* gammaStrut = new G4Box("gammaStrut_S",228.6/2.0,12.7/2.0,12.7/2.0);
      
        //make holes around
        shift = G4ThreeVector(0.0,0.0,0.0);
        gammaStrut = new G4SubtractionSolid("gammaStrut_S",gammaStrut,gammastrut_center_hole,NULL,shift);

	shift.setX(-((228.6/2.0)-4.763));
        gammaStrut = new G4SubtractionSolid("gammaStrut_S",gammaStrut,gammastrut_through_hole,NULL,shift);

	shift.setX(((228.6/2.0)-4.763));
        gammaStrut = new G4SubtractionSolid("gammaStrut_S",gammaStrut,gammastrut_through_hole,NULL,shift);

	//make notches
	shift.setX(-((228.6/2.0)+(12.7/2.0)-9.525));
	shift.setZ(-(12.7/2.0));
        gammaStrut = new G4SubtractionSolid("gammaStrut_S",gammaStrut,gammastrut_cutout,NULL,shift);
	shift.setX(+((228.6/2.0)+(12.7/2.0)-9.525));
	shift.setZ(-(12.7/2.0));
        gammaStrut = new G4SubtractionSolid("gammaStrut_S",gammaStrut,gammastrut_cutout,NULL,shift);
        /****************************End Gamma Strut Construction******************************************/

        /****************************Set up Gamma Bar Construction*****************************************/
        //gamma bar subtractions 
        G4Tubs* source_hole = new G4Tubs("source_hole",0.0,4.0/2.0,7.0/2.0,0.0,2*pi);

        //gamma bar
        G4VSolid* gammaBar = new G4Tubs("gammaBar_S",0.0,12.7/2.0,92.075/2.0,0.0,2*pi);

        //make holes around
        shift = G4ThreeVector(0.0,0.0,0.0);

	shift.setX(-((92.075/2.0)+(7.0/2.0)-(7.0-2.0))); //about 2mm sticking out
        gammaBar = new G4SubtractionSolid("gammaBar_S",gammaBar,source_hole,NULL,shift);

        /****************************End Gamma Bar Construction********************************************/

        /****************************Set up Be Plate Strut Construction************************************/
        //Be strut subtractions
        G4Tubs* bestrut_through_hole = new G4Tubs("bestrut_through_hole",0.0,4.318/2.0,25.4,0.0,2*pi);
        G4Tubs* bestrut_plate_hole = new G4Tubs("bestrut_plate_hole",0.0,4.32/2.0,25.4,0.0,2*pi);
	G4Box* bestrut_cutout_side = new G4Box("bestrut_cutout_side",12.7/2.0,25.4,12.7/2.0);
	G4Box* bestrut_cutout_mid = new G4Box("bestrut_cutout_mid",209.55/2.0,25.4,25.4);

        //Be strut
	G4VSolid* beStrut = new G4Box("beStrut_S",228.6/2.0,12.7/2.0,12.7/2.0);
      
        //make holes around
        shift = G4ThreeVector(0.0,0.0,0.0);

	shift.setX(-((228.6/2.0)-4.763));
        beStrut = new G4SubtractionSolid("beStrut_S",beStrut,bestrut_through_hole,NULL,shift);

	shift.setX(((228.6/2.0)-4.763));
        beStrut = new G4SubtractionSolid("beStrut_S",beStrut,bestrut_through_hole,NULL,shift);

	//have to make a rotation for the plate holes
        G4RotationMatrix *beplateholerot = new G4RotationMatrix;
	beplateholerot->rotateX(M_PI/2.0*rad);

	shift.setX(((209.55/2.0)-19.05));
        beStrut = new G4SubtractionSolid("beStrut_S",beStrut,bestrut_plate_hole,beplateholerot,shift);

	shift.setX(-((209.55/2.0)-19.05));
        beStrut = new G4SubtractionSolid("beStrut_S",beStrut,bestrut_plate_hole,beplateholerot,shift);

	//make notches
	shift.setX(-((228.6/2.0)+(12.7/2.0)-9.525));
	shift.setZ(-(12.7/2.0));
        beStrut = new G4SubtractionSolid("beStrut_S",beStrut,bestrut_cutout_side,NULL,shift);
	
	shift.setX(+((228.6/2.0)+(12.7/2.0)-9.525));
	shift.setZ(-(12.7/2.0));
        beStrut = new G4SubtractionSolid("beStrut_S",beStrut,bestrut_cutout_side,NULL,shift);

	shift.setX(0.0);
	shift.setY(-(25.4));
	shift.setZ(0.0);
        beStrut = new G4SubtractionSolid("beStrut_S",beStrut,bestrut_cutout_mid,NULL,shift);
        
        /****************************End Be Plate Strut Construction***************************************/

        /****************************Set up Be Source Plate Construction***********************************/
        //check some stuff
	if(BeR>=(1.5*inch - 5))
          BeR=(1.5*inch - 5.0);

	if(Bet>=(0.25*inch))
	  Bet = 0.25*inch-1.0;

        //bs plate subtractions
	//remember BeR and Bet are Beryllium radius and thickness
        G4Tubs* large_bsp_hole = new G4Tubs("large_bsp_hole",0.0,BeR-2.0,25.4,0.0,2*pi);
        G4Tubs* counter_bsp_hole = new G4Tubs("counter_bsp_hole",0.0,BeR,25.4,0.0,2*pi);
        G4Tubs* bsp_small_hole = new G4Tubs("bsp_small_hole",0.0,4.318/2.0,25.4,0.0,2*pi);
        G4Tubs* bsp_8_32tap_hole = new G4Tubs("bsp_8-32tap_hole",0.0,3.45/2.0,19.05,0.0,2*pi);

        //bs plate
	G4VSolid* bsPlate = new G4Box("bsPlate_S",207.963/2.0,211.1375/2.0,6.35/2.0);
      
        //make holes around
        shift = G4ThreeVector(0.0,0.0,0.0);

	shift.setY(-((211.1375/2.0)-103.187));
        bsPlate = new G4SubtractionSolid("bsPlate_S",bsPlate,large_bsp_hole,NULL,shift);

	shift.setY(-((211.1375/2.0)-103.187));
	shift.setZ(25.4+(6.35/2.0)-Bet); //Bet deep counterbore to hold 2mm thick wafer
        bsPlate = new G4SubtractionSolid("bsPlate_S",bsPlate,counter_bsp_hole,NULL,shift);

	shift.setY(-((211.1375/2.0)-103.187));
	shift.setZ(0.0);
	shift.setX(38.1);
        bsPlate = new G4SubtractionSolid("bsPlate_S",bsPlate,bsp_8_32tap_hole,NULL,shift);

	shift.setX(-38.1);
        bsPlate = new G4SubtractionSolid("bsPlate_S",bsPlate,bsp_8_32tap_hole,NULL,shift);

	shift.setX(0.0);
	shift.setY(-((211.1375/2.0)-103.187)+38.1);
        bsPlate = new G4SubtractionSolid("bsPlate_S",bsPlate,bsp_8_32tap_hole,NULL,shift);

	shift.setY(-((211.1375/2.0)-103.187)-38.1);
        bsPlate = new G4SubtractionSolid("bsPlate_S",bsPlate,bsp_8_32tap_hole,NULL,shift);

	shift.setY(((211.1375/2.0)-6.35));
	shift.setX(-85.72);
        bsPlate = new G4SubtractionSolid("bsPlate_S",bsPlate,bsp_small_hole,NULL,shift);

	shift.setX(+85.72);
        bsPlate = new G4SubtractionSolid("bsPlate_S",bsPlate,bsp_small_hole,NULL,shift);
        /****************************End Be Source Plate Plate Construction********************************/

        /****************************Set up Source Capsule Construction************************************/
        //source head 
        G4Tubs* sHead = new G4Tubs("sHead_S",0.0,6.4/2.0,10.6/2.0,0.0,2*pi);
        //source shaft 
        G4Tubs* sShaft = new G4Tubs("sShaft_S",0.0,4.0/2.0,7.0/2.0,0.0,2*pi);

        //join pieces
        shift = G4ThreeVector(0.0,0.0,0.0);
	shift.setZ(+((10.6/2.0)+(7.0/2.0)));
        G4VSolid* sCapsule = new G4UnionSolid("bsPlate_S",sHead,sShaft,NULL,shift);
        /****************************End Source Capsule Construction***************************************/

        /****************************Set up Be Source Wafer Construction***********************************/
        
	//be source wafer
        G4VSolid* beWafer = new G4Tubs("beWafer_S",0.0,BeR,Bet/2.0,0.0,2*pi);
        /****************************End Be Source Wafer Construction**************************************/

        /****************************Set up Lead Bricks Construction***************************************/
        //lead brick defs
	G4double brick_cvr_thk=0.1; //100 um aluminum foil

        //bricks
	G4VSolid* innerBrick = new G4Box("innerBrick_S",101.6/2.0,50.8/2.0,203.2/2.0);
	G4VSolid* outerBrick = new G4Box("outerBrick_S",(101.6/2.0)+brick_cvr_thk,(50.8/2.0)+brick_cvr_thk,(203.2/2.0)+brick_cvr_thk);
        /****************************End Lead Bricks Construction******************************************/

        /****************************Set up Bridge Construction********************************************/

	//make a cross section
        std::vector<G4TwoVector> outerVerts;
	if(doFancy){
	  //make a version that is very detailed
          //read bridge cross section file
	  std::ifstream thefile(file8020.c_str(),std::ios::in);
	  G4String firstline;
	  getline(thefile,firstline);
	  G4int count=0;
	  while(!thefile.eof()){
            G4double x,y;
	    thefile >> x >> y;
	    //convert to mm
	    x*=inch;
	    y*=inch;
	    //woops, I didn't use the right coordinates for the data theif, scale by 2
	    x/=2.0;
	    y/=2.0;
	    //if(count%10==0)
	    if(std::abs(x)<100.0){ //some reason last point is anomalous
	      outerVerts.push_back(G4TwoVector(x,y));
	      G4cout << "X: " << x << " mm; Y: " << y << " mm" << G4endl;
	    }
	    count++;
          }
	}
	else{
          outerVerts.push_back(G4TwoVector(-1.5*inch/2.0,+1.5*inch/2.0));
          outerVerts.push_back(G4TwoVector(-1.5*inch/2.0,-1.5*inch/2.0));
          outerVerts.push_back(G4TwoVector(+1.5*inch/2.0,-1.5*inch/2.0));
          outerVerts.push_back(G4TwoVector(+1.5*inch/2.0,+1.5*inch/2.0));
        }

	//hole
        G4Tubs* strut_center_hole = new G4Tubs("strut_center_hole",0.0,0.25*inch/2.0,35.0*inch/2.0,0.0,2*pi);

	//corner plate cut
	G4VSolid *corner_cut = new G4Box("corner_cut",(sqrt(2)*4.5*inch/2.0),(0.5*inch/2.0),(sqrt(2)*4.5*inch/2.0)); //contrived to have size of diagonal of plate

        //bridge components
        G4VSolid *long_bridge_strut = new G4ExtrudedSolid("long_bridge_strut",outerVerts,(32.5*inch/2.0),0.0,1.0,0.0,1.0);
        G4VSolid *short_bridge_strut = new G4ExtrudedSolid("short_bridge_strut",outerVerts,(6.5*inch/2.0),0.0,1.0,0.0,1.0);
	G4VSolid *corner_plate = new G4Box("corner_plate",(4.5*inch/2.0),(0.25*inch/2.0),(4.5*inch/2.0));

	//have to make a rotation for the corner plate cut 
        G4RotationMatrix *cornercutrot = new G4RotationMatrix;
	cornercutrot->rotateY(-M_PI/4.0*rad);

	//put a hole in the strut(s)
        shift = G4ThreeVector(0.0,0.0,0.0);
        long_bridge_strut = new G4SubtractionSolid("longBridgeStrut_S",long_bridge_strut,strut_center_hole,NULL,shift);
        short_bridge_strut = new G4SubtractionSolid("shortBridgeStrut_S",short_bridge_strut,strut_center_hole,NULL,shift);
	G4double dstar = 4.5*inch*pow(cos(M_PI/4.0),2) + 1.5*inch*pow(sin(M_PI/4.0),2);
	shift.setX(dstar);
	shift.setZ(dstar);
	corner_plate = new G4SubtractionSolid("corner_plate",corner_plate,corner_cut,cornercutrot,shift);
        /****************************End Bridge Construction*********************************************/

	if(bcpOn){
          //need a rotation placing the bridge coupling plates, about x 90 deg in the positive sense 
          G4RotationMatrix *bcplaterot = new G4RotationMatrix;
          bcplaterot->rotateX((M_PI/2.)*rad);

          //compose the rotation
	  //bcplaterot->transform(*overallrot);
          
	  //place the bridge coupling 
	  G4LogicalVolume* logicalBCPlate = new G4LogicalVolume(bcPlate,Aluminium,"bcPlate_L",0,0,0);
	  //G4VPhysicalVolume* bcPlateSrcWorld = new G4PVPlacement(&bcplaterot->transform(*overallrot), 
	  //G4VPhysicalVolume* bcPlateSrcWorld = new G4PVPlacement(&(*bcplaterot*(*overallrot)), 
	  G4VPhysicalVolume* bcPlateSrcWorld = new G4PVPlacement(bcplaterot, 
	  							//disp+G4ThreeVector(0,((203.2/2.0)+(6.35/2.0)),((241.3/2.0)-19.05+(6.35/2.0))).transform(*overallrot),
	  							disp+G4ThreeVector(0,((203.2/2.0)+(6.35/2.0)),((241.3/2.0)-19.05+(6.35/2.0))),
	  							"bcPlate_P",
	  							logicalBCPlate,
	  							world,
	  							false,
	  							0);

	  // Visualization attributes
	  G4VisAttributes* VisAttBCPlate = new G4VisAttributes(G4Colour(0.75,0.75,0.75));
	  VisAttBCPlate->SetForceWireframe(false);  //I want a Wireframe of the me
	  logicalBCPlate->SetVisAttributes(VisAttBCPlate);  
	  // Make Invisible
	  if(!bcpVis)
	    logicalBCPlate->SetVisAttributes(G4VisAttributes::Invisible);
	}

	if(bpOn){
          //place the bottom plate
	  G4LogicalVolume* logicalBotPlate = new G4LogicalVolume(botPlate,Aluminium,"botPlate_L",0,0,0);
	  G4VPhysicalVolume* botPlateSrcWorld = new G4PVPlacement(overallrot, 
	  							disp,
	  							"botPlate_P",
	  							logicalBotPlate,
	  							world,
	  							false,
	  							0);

	  // Visualization attributes
	  G4VisAttributes* VisAttBotPlate = new G4VisAttributes(G4Colour(192/255.0,192/255.0,192/255.0));
	  VisAttBotPlate->SetForceWireframe(false);  //I want a Wireframe of the me
	  logicalBotPlate->SetVisAttributes(VisAttBotPlate);  
	  // Make Invisible
	  if(!bpVis)
	    logicalBotPlate->SetVisAttributes(G4VisAttributes::Invisible);
	}

	if(spOn){
          //need a rotation placing the side plates, about y 90 deg in the negative sense 
          G4RotationMatrix *sideplaterot;
	  row1 = G4ThreeVector(0,0,1);
	  row2 = G4ThreeVector(0,1.0,0);
	  row3 = G4ThreeVector(-1,0,0);
	  sideplaterot = new G4RotationMatrix(row1,row2,row3);
	  sideplaterot->setRows(row1,row2,row3);

	  //compose with overall rotation
	  //sideplaterot->transform(*overallrot);
	  //overallrot->transform(*sideplaterot);

          //place the side plates	
	  G4LogicalVolume* logicalSidePlate = new G4LogicalVolume(sidePlate,Aluminium,"sidePlate_L",0,0,0);
	  G4VPhysicalVolume* sideLPlateSrcWorld = new G4PVPlacement(sideplaterot, 
	  							disp+G4ThreeVector((corner_x+(9.525/2.0)),0,(203.2/2.0)+(6.35/2.0)), 
	  							"sideLPlate_P",
	  							logicalSidePlate,
	  							world,
	  							false,
	  							0);

	  G4VPhysicalVolume* sideRPlateSrcWorld = new G4PVPlacement(sideplaterot, 
	  							disp+G4ThreeVector(-(corner_x+(9.525/2.0)),0,(203.2/2.0)+(6.35/2.0)),
	  							"sideRPlate_P",
	  							logicalSidePlate,
	  							world,
	  							false,
	  							0);

	  // Visualization attributes
	  G4VisAttributes* VisAttSidePlate = new G4VisAttributes(G4Colour(192/255.0,192/255.0,192/255.0));
	  VisAttSidePlate->SetForceWireframe(false);  //I want a Wireframe of the me
	  logicalSidePlate->SetVisAttributes(VisAttSidePlate);  
	  // Make Invisible
	  if(!spVis)
	    logicalSidePlate->SetVisAttributes(G4VisAttributes::Invisible);
	}
        
	if(gsOn){
          //place the gamma strut 
	  //FIXME this kinda depends on which slot we want things in -- that should be an input assume second to bottom
	  G4LogicalVolume* logicalGStrut = new G4LogicalVolume(gammaStrut,Aluminium,"gammaStrut_L",0,0,0);
	  G4VPhysicalVolume* gammaStrutWorld = new G4PVPlacement(thisrot, 
	  							disp+G4ThreeVector(0,(corner_y+(6*12.7)+(6.35/2.0)),((6.35/2.0)+(203.2))),
	  							"gammaStrut_P",
	  							logicalGStrut,
	  							world,
	  							false,
	  							0);

	  // Visualization attributes
	  G4VisAttributes* VisAttGStrut = new G4VisAttributes(G4Colour(0.75,0.75,0.75));
	  VisAttGStrut->SetForceWireframe(false);  //I want a Wireframe of the me
	  logicalGStrut->SetVisAttributes(VisAttGStrut);  
	  // Make Invisible
	  if(!gsVis)
	    logicalGStrut->SetVisAttributes(G4VisAttributes::Invisible);
	}
	 
	if(gbOn){
          //place the gamma bar 
	  //FIXME this kinda depends on which slot we want things in -- that should be an input assume second to bottom
	  G4LogicalVolume* logicalGBar = new G4LogicalVolume(gammaBar,Aluminium,"gammaBar_L",0,0,0);
	  G4VPhysicalVolume* gammaBarWorld = new G4PVPlacement(thisrot, 
	  							disp+G4ThreeVector(0,(corner_y+(6*12.7)+(6.35/2.0)),((6.35/2.0)+(203.2)-(92.075/2.0)-6.35)),
	  							"gammaBar_P",
	  							logicalGBar,
	  							world,
	  							false,
	  							0);

	  // Visualization attributes
	  G4VisAttributes* VisAttGBar = new G4VisAttributes(G4Colour(192/255.0,192/255.0,192/255.0));
	  VisAttGBar->SetForceWireframe(false);  //I want a Wireframe of the me
	  logicalGBar->SetVisAttributes(VisAttGBar);  
	  // Make Invisible
	  if(!gbVis)
	    logicalGBar->SetVisAttributes(G4VisAttributes::Invisible);
        }

	if(besOn){
          //place the Be plate strut 
	  //FIXME this kinda depends on which slot we want things in -- that should be an input assume bottom
	  G4LogicalVolume* logicalBeStrut = new G4LogicalVolume(beStrut,Aluminium,"beStrut_L",0,0,0);
	  G4VPhysicalVolume* beStrutWorld = new G4PVPlacement(thisrot, 
	  							disp+G4ThreeVector(0,(corner_y+(7*12.7)+(6.35/2.0)),((6.35/2.0)+(203.2))),
	  							"beStrut_P",
	  							logicalBeStrut,
	  							world,
	  							false,
	  							0);

	  // Visualization attributes
	  G4VisAttributes* VisAttBeStrut = new G4VisAttributes(G4Colour(192/255.0,192/255.0,192/255.0));
	  VisAttBeStrut->SetForceWireframe(false);  //I want a Wireframe of the me
	  logicalBeStrut->SetVisAttributes(VisAttBeStrut);  
	  // Make Invisible
	  if(!besVis)
	    logicalBeStrut->SetVisAttributes(G4VisAttributes::Invisible);
	}

	if(bespOn){
          //need a rotation placing the Be source plates, about x 90 deg in the positive sense 
          G4RotationMatrix *bsplaterot = new G4RotationMatrix;
          bsplaterot->rotateX(-M_PI/2.*rad);
	  
          //place the Be plate 
	  //FIXME this kinda depends on which slot we want things in -- that should be an input assume bottom
	  G4LogicalVolume* logicalBSPlate = new G4LogicalVolume(bsPlate,Aluminium,"bsPlate_L",0,0,0);
	  G4VPhysicalVolume* bsPlateWorld = new G4PVPlacement(bsplaterot, 
	  							disp+G4ThreeVector(0,(corner_y+(7*12.7)+(6.35/2.0)),((6.35/2.0)+(203.2/2.0)+2.381)),
	  							"bsPlate_P",
	  							logicalBSPlate,
	  							world,
	  							false,
	  							0);

	  // Visualization attributes
	  G4VisAttributes* VisAttBSPlate = new G4VisAttributes(G4Colour(192/255.0,192/255.0,192/255.0));
	  VisAttBSPlate->SetForceWireframe(false);  //I want a Wireframe of the me
	  logicalBSPlate->SetVisAttributes(VisAttBSPlate);  
	  // Make Invisible
	  if(!bespVis)
	    logicalBSPlate->SetVisAttributes(G4VisAttributes::Invisible);
	}
	
	if(scOn){
          //place the source capsule 
	  G4LogicalVolume* logicalSCapsule = new G4LogicalVolume(sCapsule,stainlessSteel,"sCapsule_L",0,0,0);
	  G4VPhysicalVolume* sCapsuleWorld = new G4PVPlacement(thisrot, 
	  							disp+G4ThreeVector(0,(corner_y+(6*12.7)+(6.35/2.0)),((6.35/2.0)+(203.2)-6.35-92.075-((10.6/2.0)+2.0))), //2mm stick out of gamma bar
	  							"sCapsule_P",
	  							logicalSCapsule,
	  							world,
	  							false,
	  							0);

	  // Visualization attributes
	  G4VisAttributes* VisAttSCapsule = new G4VisAttributes(G4Colour(128/255.0,128/255.0,128/255.0));
	  VisAttSCapsule->SetForceWireframe(false);  //I want a Wireframe of the me
	  logicalSCapsule->SetVisAttributes(VisAttSCapsule);  
	  // Make Invisible
	  if(!scVis)
	    logicalSCapsule->SetVisAttributes(G4VisAttributes::Invisible);
	}

	if(bswOn){
          //need a rotation placing the Be source plates, about x 90 deg in the positive sense 
          G4RotationMatrix *bewaferrot = new G4RotationMatrix;
          bewaferrot->rotateX(-M_PI/2.*rad);

          //place the source capsule 
	  G4LogicalVolume* logicalBeWafer;
	  if(UseBePure)
	    logicalBeWafer = new G4LogicalVolume(beWafer,BePure,"beWafer_L",0,0,0);
	  else
	    logicalBeWafer = new G4LogicalVolume(beWafer,BeO,"beWafer_L",0,0,0);
	  G4VPhysicalVolume* beWaferWorld = new G4PVPlacement(bewaferrot, 
	  							disp+G4ThreeVector(0,(corner_y+(7*12.7)+(Bet/2.0)),(6.35/2.0)+(203.2/2.0)), 
	  							"beWafer_P",
	  							logicalBeWafer,
	  							world,
	  							false,
	  							0);

	  // Visualization attributes
	  G4VisAttributes* VisAttBeWafer = new G4VisAttributes(G4Colour(1.0,0,0));
	  VisAttBeWafer->SetForceWireframe(false);  //I want a Wireframe of the me
	  logicalBeWafer->SetVisAttributes(VisAttBeWafer);  
	  // Make Invisible
	  if(!bswVis)
	    logicalBeWafer->SetVisAttributes(G4VisAttributes::Invisible);
	}

	if(lbrickOn){
          //place the brick(s) 
	  G4LogicalVolume* logical10Brick = new G4LogicalVolume(outerBrick,Aluminium,"10Brick_L",0,0,0);
	  G4VPhysicalVolume* brick10World = new G4PVPlacement(thisrot, 
	  							disp+G4ThreeVector((101.6/2.0)+brick_cvr_thk,(101.6)-(50.6/2.0)-brick_cvr_thk,(203.2/2.0)+brick_cvr_thk), 
	  							"10Brick_P",
	  							logical10Brick,
	  							world,
	  							false,
	  							0);
	  G4LogicalVolume* logical10InnerBrick = new G4LogicalVolume(innerBrick,Pb,"10InnerBrick_L",0,0,0);
	  G4VPhysicalVolume* brick10InnerWorld = new G4PVPlacement(0, 
	  							G4ThreeVector(0,0,0), 
	  							"10InnerBrick_P",
	  							logical10InnerBrick,
	  							brick10World,
	  							false,
	  							0);

	  G4LogicalVolume* logical11Brick = new G4LogicalVolume(outerBrick,Aluminium,"11Brick_L",0,0,0);
	  G4VPhysicalVolume* brick11World = new G4PVPlacement(thisrot, 
	  							disp+G4ThreeVector(-((101.6/2.0)+brick_cvr_thk),(101.6)-(50.6/2.0)-brick_cvr_thk,(203.2/2.0)+brick_cvr_thk), 
	  							"11Brick_P",
	  							logical11Brick,
	  							world,
	  							false,
	  							0);
	  G4LogicalVolume* logical11InnerBrick = new G4LogicalVolume(innerBrick,Pb,"11InnerBrick_L",0,0,0);
	  G4VPhysicalVolume* brick11InnerWorld = new G4PVPlacement(0, 
	  							G4ThreeVector(0,0,0), 
	  							"11InnerBrick_P",
	  							logical11InnerBrick,
	  							brick11World,
	  							false,
	  							0);

          //next brick layer needs a rotation, about z 90 deg in the positive sense 
          G4RotationMatrix *brickrot = new G4RotationMatrix;
          brickrot->rotateY(-M_PI/2.*rad);

          //place the brick(s) upper level
	  G4LogicalVolume* logical00Brick = new G4LogicalVolume(outerBrick,Aluminium,"00Brick_L",0,0,0);
	  G4VPhysicalVolume* brick00World = new G4PVPlacement(brickrot, 
	  							disp+G4ThreeVector(0,(203.2/2.0)-(50.6/2.0)-brick_cvr_thk-(50.6)-2*brick_cvr_thk,+(6.35/2.0)+(101.6/2.0)+brick_cvr_thk+101.6+2*brick_cvr_thk), 
	  							"00Brick_P",
	  							logical00Brick,
	  							world,
	  							false,
	  							0);
	  G4LogicalVolume* logical00InnerBrick = new G4LogicalVolume(innerBrick,Pb,"00InnerBrick_L",0,0,0);
	  G4VPhysicalVolume* brick00InnerWorld = new G4PVPlacement(0, 
	  							G4ThreeVector(0,0,0), 
	  							"00InnerBrick_P",
	  							logical00InnerBrick,
	  							brick00World,
	  							false,
	  							0);

	  G4LogicalVolume* logical01Brick = new G4LogicalVolume(outerBrick,Aluminium,"01Brick_L",0,0,0);
	  G4VPhysicalVolume* brick01World = new G4PVPlacement(brickrot, 
	  							disp+G4ThreeVector(0,(203.2/2.0)-(50.6/2.0)-brick_cvr_thk-(50.6)-2*brick_cvr_thk,+(6.35/2.0)+(101.6/2.0)+brick_cvr_thk), 
	  							"01Brick_P",
	  							logical01Brick,
	  							world,
	  							false,
	  							0);
	  G4LogicalVolume* logical01InnerBrick = new G4LogicalVolume(innerBrick,Pb,"01InnerBrick_L",0,0,0);
	  G4VPhysicalVolume* brick01InnerWorld = new G4PVPlacement(0, 
	  							G4ThreeVector(0,0,0), 
	  							"01InnerBrick_P",
	  							logical01InnerBrick,
	  							brick01World,
	  							false,
	  							0);
	  // Visualization attributes
	  G4VisAttributes* VisAttBricks = new G4VisAttributes(G4Colour(255/255.0,153/255.0,51/255.0));
	  VisAttBricks->SetForceWireframe(false);  //I want a Wireframe of the me
	  logical10Brick->SetVisAttributes(VisAttBricks);  
	  logical11Brick->SetVisAttributes(VisAttBricks);  
	  logical00Brick->SetVisAttributes(VisAttBricks);  
	  logical01Brick->SetVisAttributes(VisAttBricks);  
	  // Make Invisible
	  if(!lbrickVis){
	    logical10Brick->SetVisAttributes(G4VisAttributes::Invisible);
	    logical11Brick->SetVisAttributes(G4VisAttributes::Invisible);
	    logical00Brick->SetVisAttributes(G4VisAttributes::Invisible);
	    logical01Brick->SetVisAttributes(G4VisAttributes::Invisible);
	    logical10InnerBrick->SetVisAttributes(G4VisAttributes::Invisible);
	    logical11InnerBrick->SetVisAttributes(G4VisAttributes::Invisible);
	    logical00InnerBrick->SetVisAttributes(G4VisAttributes::Invisible);
	    logical01InnerBrick->SetVisAttributes(G4VisAttributes::Invisible);
	  }
	}
	
	if(bstructOn){
          //need a rotation placing the bridge struts, about x 90 deg in the positive sense 
          G4RotationMatrix *bridgestrutrot = new G4RotationMatrix;
          bridgestrutrot->rotateY(M_PI/2.*rad);

          //place the struts 
	  G4LogicalVolume* logicalLongBridgeStrut = new G4LogicalVolume(long_bridge_strut,Aluminium,"longBridgeStrut_L",0,0,0);
	  G4VPhysicalVolume* longBridgeStrutBack = new G4PVPlacement(bridgestrutrot, 
	  							disp+G4ThreeVector(0,(1.5*inch/2.0)+(4.0*inch)+(0.25*inch),(1.5*inch/2.0)+(0.25*inch/2.0)+(4.0*inch)-(9.5*inch/2.0)),
	  							"longBridgeStrutBack_P",
	  							logicalLongBridgeStrut,
	  							world,
	  							false,
	  							0);
	  G4VPhysicalVolume* longBridgeStrutFront = new G4PVPlacement(bridgestrutrot, 
	  							disp+G4ThreeVector(0,(1.5*inch/2.0)+(4.0*inch)+(0.25*inch),-(1.5*inch/2.0)+(0.25*inch/2.0)+(4.0*inch)+(9.5*inch/2.0)),
	  							"longBridgeStrutFront_P",
	  							logicalLongBridgeStrut,
	  							world,
	  							false,
	  							0);


	  //place the short struts
	  G4LogicalVolume* logicalShortBridgeStrut = new G4LogicalVolume(short_bridge_strut,Aluminium,"shortBridgeStrut_L",0,0,0);
	  G4VPhysicalVolume* shortBridgeStrutLeft = new G4PVPlacement(thisrot, 
	  							disp+G4ThreeVector(-((1.5*inch/2.0)-(32.5*inch/2)),(1.5*inch/2.0)+(4.0*inch)+(0.25*inch),(0.25*inch/2.0)+(4.0*inch)),
	  							"shortBridgeStrutLeft_P",
	  							logicalShortBridgeStrut,
	  							world,
	  							false,
	  							0);
	  G4VPhysicalVolume* shortBridgeStrutRight = new G4PVPlacement(thisrot, 
	  							disp+G4ThreeVector(((1.5*inch/2.0)-(32.5*inch/2.0)),(1.5*inch/2.0)+(4.0*inch)+(0.25*inch),(0.25*inch/2.0)+(4.0*inch)),
	  							"shortBridgeStrutRight_P",
	  							logicalShortBridgeStrut,
	  							world,
	  							false,
	  							0);
	  G4VPhysicalVolume* shortBridgeStrutLeftClose = new G4PVPlacement(thisrot, 
	  							disp+G4ThreeVector(-((1.5*inch/2.0)-(11.5*inch/2.0)),(1.5*inch/2.0)+(4.0*inch)+(0.25*inch),(0.25*inch/2.0)+(4.0*inch)),
	  							"shortBridgeStrutRight_P",
	  							logicalShortBridgeStrut,
	  							world,
	  							false,
	  							0);
	  G4VPhysicalVolume* shortBridgeStrutRightClose = new G4PVPlacement(thisrot, 
	  							disp+G4ThreeVector(((1.5*inch/2.0)-(11.5*inch/2.0)),(1.5*inch/2.0)+(4.0*inch)+(0.25*inch),(0.25*inch/2.0)+(4.0*inch)),
	  							"shortBridgeStrutRight_P",
	  							logicalShortBridgeStrut,
	  							world,
	  							false,
	  							0);

	  //place the corner plates 
	  G4LogicalVolume* logicalBridgeCornerPlate = new G4LogicalVolume(corner_plate,Aluminium,"bridgeCornerPlate_L",0,0,0);
	  G4VPhysicalVolume* bridgeCornerPlateRightBack = new G4PVPlacement(thisrot, 
	  							disp+G4ThreeVector(((4.5*inch/2.0)-(32.5*inch/2)),+(0.25*inch/2.0)+(4.0*inch),(4.5*inch/2.0)+(0.25*inch/2.0)+(4*inch)-(9.5*inch/2.0)),
	  							"bridgeCornerPlateRightBack_P",
	  							logicalBridgeCornerPlate,
	  							world,
	  							false,
	  							0);

          //need a rotation z by 180 deg
          G4RotationMatrix *leftbackrot = new G4RotationMatrix;
          leftbackrot->rotateZ(M_PI*rad);

	  G4VPhysicalVolume* bridgeCornerPlateLeftBack = new G4PVPlacement(leftbackrot, 
	  							disp+G4ThreeVector(-((4.5*inch/2.0)-(32.5*inch/2)),+(0.25*inch/2.0)+(4.0*inch),(4.5*inch/2.0)+(0.25*inch/2.0)+(4*inch)-(9.5*inch/2.0)),
	  							"bridgeCornerPlateLeftBack_P",
	  							logicalBridgeCornerPlate,
	  							world,
	  							false,
	  							0);

          //need a rotation x by 180 deg, then rotate z by 180 deg
          G4RotationMatrix *leftfrontrot = new G4RotationMatrix;
          leftfrontrot->rotateX(M_PI*rad);
          leftfrontrot->rotateZ(M_PI*rad);

	  G4VPhysicalVolume* bridgeCornerPlateLeftFront = new G4PVPlacement(leftfrontrot, 
	  							disp+G4ThreeVector(-((4.5*inch/2.0)-(32.5*inch/2)),+(0.25*inch/2.0)+(4.0*inch),-(4.5*inch/2.0)+(0.25*inch/2.0)+(4*inch)+(9.5*inch/2.0)),
	  							"bridgeCornerPlateLeftFront_P",
	  							logicalBridgeCornerPlate,
	  							world,
	  							false,
	  							0);

          //need a rotation x by 180 deg
          G4RotationMatrix *rightfrontrot = new G4RotationMatrix;
          rightfrontrot->rotateX(M_PI*rad);

	  G4VPhysicalVolume* bridgeCornerPlateRightFront = new G4PVPlacement(rightfrontrot, 
	  							disp+G4ThreeVector(+((4.5*inch/2.0)-(32.5*inch/2)),+(0.25*inch/2.0)+(4.0*inch),-(4.5*inch/2.0)+(0.25*inch/2.0)+(4*inch)+(9.5*inch/2.0)),
	  							"bridgeCornerPlateRightFront_P",
	  							logicalBridgeCornerPlate,
	  							world,
	  							false,
	  							0);
	  // Visualization attributes
	  G4VisAttributes* VisAttBridge = new G4VisAttributes(G4Colour(192/255.0,192/255.0,192/255.0));
	  G4VisAttributes* VisAttCrossStrut = new G4VisAttributes(G4Colour(1.0,0,0));
	  VisAttBridge->SetForceWireframe(false);  //I want a Wireframe of the me
	  VisAttCrossStrut->SetForceWireframe(false);  //I want a Wireframe of the me
	  logicalLongBridgeStrut->SetVisAttributes(VisAttBridge);  
	  logicalShortBridgeStrut->SetVisAttributes(VisAttBridge);  
	  logicalBridgeCornerPlate->SetVisAttributes(VisAttBridge);  
	  // Make Invisible
	  if(!bstructVis){
	    logicalLongBridgeStrut->SetVisAttributes(G4VisAttributes::Invisible);
	    logicalShortBridgeStrut->SetVisAttributes(G4VisAttributes::Invisible);
	    logicalBridgeCornerPlate->SetVisAttributes(G4VisAttributes::Invisible);
	  }
        }

        return;
}
void NeutReflect_DetectorConstruction::ConstructPhotoNSourceBoxNested(G4VPhysicalVolume *world,G4double xpos,
G4double ypos,
G4double zpos,
G4int bshelf,
G4int gshelf)
{
        //only use rotation if parent volume is null (highest)
        G4RotationMatrix *thisrot;
        G4cout << world->GetName() << G4endl;
        if(world->GetName()=="world_P"){
          thisrot=xrot;
        }
        else{
          thisrot=nullrot;
        }

        //displacement vector
        G4ThreeVector disp = G4ThreeVector(xpos,ypos,zpos).transform(*thisrot);

        //various geometry objects needed for computations
        G4ThreeVector shift,row1,row2,row3;

        /****************************Set up Bridge Coupling Plate Construction***********************************/

	//bridge coupling plate
	G4VSolid* bcPlate = new G4Box("bcPlate_S",203.2/2.0,203.2/2.0,6.35/2.0);

        /****************************End Bridge Coupling Plate Construction**************************************/


        /****************************Set up Bottom Plate Construction*******************************************/
	/*
        //bottom plate subtractions
        G4Tubs* small_hole = new G4Tubs("small_hole",0.0,4.32/2.0,10,0.0,2*pi);
        G4Tubs* large_hole = new G4Tubs("large_hole",0.0,4.32/2.0,10,0.0,2*pi);
        G4Box* shelf = new G4Box("shelf",209.55/2.0,6.35/2.0,10.0/2.0);

        //bottom plate
	G4VSolid* botPlate = new G4Box("botPlate_S",228.6/2.0,203.2/2.0,6.35/2.0);

        //make holes around
        float corner_x=-114.3,corner_y=-101.6;
        float ystep=25.4;
        float yridgestep=15.88;
        float x1=4.76,x2=223.8;
        shift = G4ThreeVector(0.0,0.0,0.0);
        for(int i=0;i<7;i++){

           //get the iteration
      	   std::ostringstream oss;
	   oss << i;
        
           //small holes on the "left"
           shift.setX(corner_x+x1);
           shift.setY(corner_y+(i+1)*ystep);
           botPlate = new G4SubtractionSolid("botPlate_S_left_"+oss.str(),botPlate,small_hole,NULL,shift);

           //small holes on the "right"
           shift.setX(corner_x+x2);
           shift.setY(corner_y+(i+1)*ystep);
           botPlate = new G4SubtractionSolid("botPlate_S_right_"+oss.str(),botPlate,small_hole,NULL,shift);
        }
        for(int i=0;i<8;i++){
        
           //get the iteration
	   std::ostringstream oss;
	   oss << i;

           //ridges
           shift.setX(0.0);
           shift.setY(corner_y+(i+1)*yridgestep);
           shift.setZ((10.0/2.0)-1.59); //1/16" depth
           botPlate = (G4VSolid*) new G4SubtractionSolid("botPlate_S_ridge_"+oss.str(),botPlate,shelf,NULL,shift);

        }
        */
        /****************************End Bottom Plate Construction*******************************************/

        /****************************Set up Side Plate Construction*******************************************/
        /*//side plate subtractions
        G4Tubs* small_side_hole = new G4Tubs("small_side_hole",0.0,3.45/2.0,19.05,0.0,2*pi);

        //side plate
	G4VSolid* sidePlate = new G4Box("sidePlate_S",203.2/2.0,203.2/2.0,6.35/2.0);

        //need a rotation for these holes
        G4RotationMatrix *holerot;
	row1 = G4ThreeVector(1.0,0,0);
	row2 = G4ThreeVector(0,0,1.0);
	row3 = G4ThreeVector(0,-1.0,0);
	holerot = new G4RotationMatrix(row1,row2,row3);
	holerot->setRows(row1,row2,row3);

        //make holes around
        float s_corner_x=-101.6,s_corner_y=-101.6;
        float s_ystep=25.4;
        shift = G4ThreeVector(0.0,0.0,0.0);
        for(int i=0;i<7;i++){
        
           //small holes on the "left"
           shift.setX(s_corner_x);
           shift.setY(s_corner_y+(i+1)*s_ystep);
           botPlate = (G4VSolid*) new G4SubtractionSolid("botPlate_S",botPlate,small_side_hole,holerot,shift);

        }
        */
        /****************************End Side Plate Construction*******************************************/
        
        //place the bottom plate
	/*G4LogicalVolume* logicalBotPlate = new G4LogicalVolume(botPlate,Aluminium,"botPlate_L",0,0,0);
	G4VPhysicalVolume* botPlateSrcWorld = new G4PVPlacement(thisrot, 
								disp,
								"botPlate_P",
								logicalBotPlate,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttBotPlate = new G4VisAttributes(G4Colour(0.75,0.75,0.75));
	VisAttBotPlate->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalBotPlate->SetVisAttributes(VisAttBotPlate);  
	// Make Invisible
	//logicalBotPlate->SetVisAttributes(G4VisAttributes::Invisible);
        */

        //place the side plates	
	/*G4LogicalVolume* logicalSidePlate = new G4LogicalVolume(sidePlate,Aluminium,"sidePlate_L",0,0,0);
	G4VPhysicalVolume* sideLPlateSrcWorld = new G4PVPlacement(thisrot, 
								disp+G4ThreeVector(0,0,2*25.4), // 2 inches up
								"sideLPlate_P",
								logicalSidePlate,
								world,
								false,
								0);

	G4VPhysicalVolume* sideRPlateSrcWorld = new G4PVPlacement(thisrot, 
								disp+G4ThreeVector(0,0,4*25.4), // 4 inches up
								"sideRPlate_P",
								logicalSidePlate,
								world,
								false,
								0);

	// Visualization attributes
	G4VisAttributes* VisAttSidePlate = new G4VisAttributes(G4Colour(0.75,0.75,0.75));
	VisAttSidePlate->SetForceWireframe(false);  //I want a Wireframe of the me
	logicalSidePlate->SetVisAttributes(VisAttSidePlate);  
	// Make Invisible
	//logicalSidePlate->SetVisAttributes(G4VisAttributes::Invisible);
        */
	
        return;
}
// ------------------Tracking Region----------------------

void NeutReflect_DetectorConstruction::ConstructTracking()
{
	
	// Sensor Region (cuts 250 eV) for any particle within the DetBox logical volume
	//  G4Region* DetectorRegion;
	
	//DetectorRegion = new G4Region(G4String("Detector"));
	//logicalDetectorBox->SetRegion(DetectorRegion);
	//DetectorRegion->AddRootLogicalVolume(logicalDetectorBox);
	
	
	
}
// ----------------Sensitive Detector ----------------------

void NeutReflect_DetectorConstruction::ConstructGenericSensitive(G4LogicalVolume* logicalGeneric,G4String name)
{
	
	
	//------------------------------------------------ 
	// Sensitive detectors setup
	//------------------------------------------------ 
	
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

        NeutReflect_StdSD* SD1 = new NeutReflect_StdSD(name,NeutReflectCollName.size()+1000);
        NeutReflectCollName[name] = NeutReflectCollName.size()+1000;
        SDman->AddNewDetector(SD1);
        logicalGeneric->SetSensitiveDetector(SD1);
    
    
}
//-------------------Routine for Updating--------------------
void NeutReflect_DetectorConstruction::UpdateGeometry()
{
	//delete DetectorRegion; // otherwise this causes segmentation faults.
	
	G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
	
}
