/*==================NeutReflect_PrimaryGeneratorAction.cc================= 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Code supporting the generators for the NeutReflect simulation
               all user-defined particle sources should be instantiated
	       in the constructor and executed in GeneratePrimaries(*).
              
              
======================================================================*/

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "Randomize.hh"
#include <math.h>

#include "NeutReflect_PrimaryGeneratorAction.hh"
#include "NeutReflect_DetectorConstruction.hh"
#include "NeutReflect_EventInfo.hh"

NeutReflect_PrimaryGeneratorAction::NeutReflect_PrimaryGeneratorAction(G4double bias, G4String dec, G4String rfile):
xnbias(bias),
decayer(dec),
rootfile(rfile)
{
  //initialize TTree pointer to NULL
  //t = NULL;

  //finalpathlength initialize
  finalpathlength = -1.0;
  
  //Get the geometry
  Ndet = (NeutReflect_DetectorConstruction*)(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  //set basic particle gun
  particleGun = new G4ParticleGun();

  //open a root file
  if(decayer=="cascade" && rootfile!="NULL"){
    f = new TFile(rootfile.c_str());
    t = (TTree*) f->Get("cascade");
    t->Print();
    t->SetBranchAddress("n",&ncascade);
    t->SetBranchAddress("time",time);
    t->SetBranchAddress("Eg",Eg);
    t->GetEntry(0);
    rootnevent=0;

    G4cout << "For zeroth entry" << G4endl;
    G4cout << "n: " << ncascade << G4endl;
    for(int i=0;i<ncascade;i++)
      G4cout << "time[" << i <<"]: " << time[i] <<  "; Eg[" << i << "]: " << Eg[0] << G4endl;
  }

  //set overall rotation to be consistent with geometry
  particleGun->SetParticlePosition(G4ThreeVector(0.0*m,0.0*m,0.0*m));
  G4ThreeVector row1 = G4ThreeVector(1,0,0);
  G4ThreeVector row2 = G4ThreeVector(0,0,1.0);
  G4ThreeVector row3 = G4ThreeVector(0,-1.0,0);
  xrot = new G4RotationMatrix(row1,row2,row3);
  xrot->setRows(row1,row2,row3);

  //set some basic kinematics
  delbe = 11.348;
  deln = 8.071;
  melec = 0.510998;
  mnbar = 931.494045;
  mbe = 9*mnbar + delbe -4*melec;
   
  
  //source types 1) gamma - isotropic gamma at x,y,z with an E spectrum consistent for sources
  //		 2) neutron -isotropic neutron at x,y,z consistent with E spectrum for gamma induced
  //                         neutrons with gammas in the 0 x - 1 y + 0 z direction x,y,z is distributed
  //                         over a disk
  //             3) neutrondome - isotropic neutron generated from a dome around 0,0,0 with some thickness
  //                              and consistent with E spectrum for gamma induced neutrons with gammas
  //                              comming from 0,0,0

  //do a different source type if decayer is "cascade"
  if(decayer=="cascade")
    sourceType="cascade";
  else
    sourceType="neutron";
  R = Ndet->GetBeR(); //radius (mm) of disk or hemisphere for neutron sources
  thk = Ndet->GetBet(); //thickness (mm) of disk or hemisphere shell for neutron sources
  //d = -1.0; //distance (mm) above or below origin for disk source (negative for below)
  //d = -3.175; //distance (mm) above or below origin for disk source (negative for below) Zeigler disk is 6.35mm thick approx
  d = Ndet->GetDiskDisplacement(); //distance (mm) above or below origin for disk source (negative for below) Zeigler disk is 6.35mm thick approx
  if(sourceType.find("gamma") != G4String::npos)
    particleGun->SetParticleDefinition(G4Gamma::Definition()); 
  else if(sourceType.find("neutron") != G4String::npos)
    particleGun->SetParticleDefinition(G4Neutron::Definition());
  else
    particleGun->SetParticleDefinition(G4Gamma::Definition()); 
}
NeutReflect_PrimaryGeneratorAction::~NeutReflect_PrimaryGeneratorAction()
{
  //if(decayer=="cascade" && rootfile=="NULL")
  //  f->Close();
  delete particleGun;
}
void NeutReflect_PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //make an ancillary EventInfo container to store some event info
  //and register it with this event
  NeutReflect_EventInfo *anEventInfo = new NeutReflect_EventInfo();
  anEvent->SetUserInformation(anEventInfo);

  //use below if you use a custom particle generator which
  //exits based on some condition (other than total number
  //of primaries thrown)
  /*if(!particleGun->GetHasExited())
   particleGun->GeneratePrimaryVertex(anEvent);
  else
  {*/
    // RunManager cannot abort the event from inside
    // UserGeneratePrimaries(), so we do a soft abort
    // to the RunManager, and abort the event ourselves.
    // The result is the same as a hard abort.
  /*     G4cout << "about to abort" << G4endl;
       G4RunManager::GetRunManager()->AbortRun(true); 
       anEvent->SetEventAborted();
  }*/

  if(sourceType=="gamma"){
    std::vector<G4double> angles;
    //set energy to come from realistic spectrum
    std::vector<G4double> genergy = GenerateGammaEnergy(decayer);
    //set the direction to an isotropic for each gamma
    for(G4int i=0;i<(G4int)genergy.size();i++){
      angles = GenerateRandomDirection();
      particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]).transform(*xrot));
      particleGun->SetParticleEnergy(genergy[i]);
      particleGun->GeneratePrimaryVertex(anEvent);
    }
    //now generate back-to-back 511's
    //this only actually happens about 0.21% of the time, else it's EC, so account for that the fraction is 0.0021 -- this doesn't get the correlation with 1.8 MeV gamma
    //correct!
    if(G4UniformRand() < 0.0021){
      angles = GenerateRandomDirection();
      particleGun->SetParticleEnergy(0.511);
      particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]).transform(*xrot));
      particleGun->GeneratePrimaryVertex(anEvent);
      particleGun->SetParticleMomentumDirection(G4ThreeVector(-angles[2]*cos(angles[0]),-angles[2]*sin(angles[0]),-angles[1]).transform(*xrot));
      particleGun->GeneratePrimaryVertex(anEvent);
    }
  }
  else if(sourceType=="neutron"){
    std::vector<G4double> angles;
    //set gamma energy to come from realistic spectrum
    std::vector<G4double> genergy = GenerateGammaEnergy(decayer);
    //set directions to isotropic for each gamma
    for(G4int i=0;i<(G4int)genergy.size();i++){
      angles = GenerateRandomDirection();
      //find out the location and thickness of the BeO struck
      std::vector<G4double> diskPath = GenerateDiskPathlength(angles);
      //std::vector<G4double> diskPath = GenerateShellPathlength(angles);
      
      //set the event information
      anEventInfo->SetBeLength(diskPath[0]);

      //generate the cross section for interaction in that BeO
      if(DidGammaNHappen(diskPath[0],genergy[i])){


        //get the precise place of interation along path (thin target approx)
        //G4cout << "throwing from distance " << diskPath[1] << " + " << diskPath[0] << G4endl;
        G4double distance = diskPath[1] + G4UniformRand()*diskPath[0];
        particleGun->SetParticlePosition(G4ThreeVector(distance*angles[2]*cos(angles[0]),distance*angles[2]*sin(angles[0]),distance*angles[1]).transform(*xrot));

        //get the angle and energy wrt the direction of the gamma
        std::vector<G4double> neutangles = GenerateRandomDirection();
        G4double delrec = 4.942;
        G4double mrec = 8*mnbar + delrec - 4*melec; //8Be is kind of like two alphas, right?
        mn = mnbar + deln; //neutron mass
        G4double Q = mbe - (mrec + mn);  //should be negative (obvi)
        //G4cout << "Q value is: " << Q << G4endl;
        G4double nenergy = GenerateGammaNEnergy(neutangles[1],genergy[i],mrec,Q);

        if(nenergy>0.0){
          //G4cout << "Neutron energy is: " << nenergy << " MeV" << G4endl;
          //create rotation to bring neutron momentup back to lab frame
          G4ThreeVector row1,row2,row3;
          row1 = G4ThreeVector(1.0,0,0.0);
          row2 = G4ThreeVector(0,angles[1],angles[2]);
          row3 = G4ThreeVector(0,-angles[2],angles[1]);
          G4RotationMatrix *neutrot1 = new G4RotationMatrix(row1,row2,row3);
          row1 = G4ThreeVector(cos(angles[0]),sin(angles[0]),0.0);
          row2 = G4ThreeVector(-sin(angles[0]),cos(angles[0]),0.0);
          row3 = G4ThreeVector(0,0,1.0);
          G4RotationMatrix *neutrot2 = new G4RotationMatrix(row1,row2,row3);

          //set up the neutron for the above angles
          particleGun->SetParticleMomentumDirection(G4ThreeVector(neutangles[2]*cos(neutangles[0]),neutangles[2]*sin(neutangles[0]),neutangles[1]).transform(*neutrot1).transform(*neutrot2).transform(*xrot));
          particleGun->SetParticleEnergy(nenergy);
          particleGun->GeneratePrimaryVertex(anEvent);
        }
        //throw the particles
      }
      else{
        //throw gamma?
        particleGun->SetParticleDefinition(G4Gamma::Definition()); 
        particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]).transform(*xrot));
        particleGun->SetParticleEnergy(genergy[i]);
        particleGun->GeneratePrimaryVertex(anEvent);
        particleGun->SetParticleDefinition(G4Neutron::Definition()); 
      }
    }
    //now generate back-to-back 511's
    //this only actually happens about 0.21% of the time, else it's EC, so account for that the fraction is 0.0021 -- this doesn't get the correlation with 1.8 MeV gamma FIXME
    //correct!
    if(G4UniformRand() < 0.0021){
      angles = GenerateRandomDirection();
      particleGun->SetParticleDefinition(G4Gamma::Definition()); 
      particleGun->SetParticleEnergy(0.511);
      particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]).transform(*xrot));
      particleGun->GeneratePrimaryVertex(anEvent);
      particleGun->SetParticleMomentumDirection(G4ThreeVector(-angles[2]*cos(angles[0]),-angles[2]*sin(angles[0]),-angles[1]).transform(*xrot));
      particleGun->GeneratePrimaryVertex(anEvent);
      particleGun->SetParticleDefinition(G4Neutron::Definition()); 
    }  
    //particle gun persists so return to appropriate value  
    particleGun->SetParticlePosition(G4ThreeVector(0.0*m,0.0*m,0.0*m));
  }
  else if(sourceType=="neutronimpl"){
    //generate a Geant4 88Y event in the source material and let Geant do the conversion in the BeO -- you must have the BeO in the sim
    //source should be at the origin in the center of what I assume to be a 5mm diameter and 3.18mm height cylinder
    G4double srcr,srcthet,srcz;
    srcthet=2*pi*G4UniformRand();
    srcz=3.18*G4UniformRand()-(3.18/2);
    //generate uniformly in r**2/2
    srcr=(2.5*2.5/2.0)*G4UniformRand();
    srcr=sqrt(2*srcr);

    //Hassan's code
    // Y-88 source
    G4int Z = 39, A = 88;
    G4double ionCharge = 0.*eplus;
    G4double excitEnergy = 0.*keV;
  
    G4ParticleDefinition* ion = 
    G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z,A,excitEnergy);
    particleGun->SetParticleDefinition(ion);
    particleGun->SetParticleCharge(ionCharge);
    particleGun->SetParticleEnergy(0.100);
    particleGun->SetParticlePosition(G4ThreeVector(srcr*cos(srcthet)*mm,srcr*sin(srcthet)*mm,srcz*mm).transform(*xrot));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(1,0,0).transform(*xrot));
    particleGun->GeneratePrimaryVertex(anEvent);

    //particle gun persists so return to appropriate value(s)  
    particleGun->SetParticlePosition(G4ThreeVector(0.0*m,0.0*m,0.0*m));
    particleGun->SetParticleDefinition(G4Gamma::Definition()); 
  }
  else if(sourceType=="cascade"){
    
    //define angle in case need it
    std::vector<G4double> angles;

    //get particle definition
    particleGun->SetParticleDefinition(G4Gamma::Definition()); 

    //set a location inside detector
    //FIXME: hard-codes the detector dimension and shape
    if(Ndet->GetDesignNo() > -9){
      G4double srcr,srcthet,srcz;
      srcthet=2*pi*G4UniformRand();
      srcz=(25.4+12.7)*G4UniformRand()-((25.4+12.7)/2);
      //generate uniformly in r**2/2
      srcr=(76.2*76.2/2.0)*G4UniformRand();
      srcr=sqrt(2*srcr);

      //set the particle position
      particleGun->SetParticlePosition(G4ThreeVector(srcr*cos(srcthet)*mm,srcr*sin(srcthet)*mm,srcz*mm).transform(*xrot));
    }
    else if(Ndet->GetDesignNo() == -9){
      G4double srcr,srcthet,srcphi;
      //generate uniform in costhet
      srcthet=-1+2*G4UniformRand();
      srcthet=acos(srcthet);
      srcphi=2*pi*G4UniformRand();
      //generate uniformly in r**3/3
      srcr=((1000.0*1000.0*1000.0)/3.0)*G4UniformRand();
      srcr=pow(3*srcr,(1/3.0));

      //set the particle position
      particleGun->SetParticlePosition(G4ThreeVector(srcr*cos(srcphi)*cos(srcthet)*mm,srcr*sin(srcphi)*cos(srcthet)*mm,srcr*sin(srcthet)*mm).transform(*xrot));
    
    }
    else if(Ndet->GetDesignNo() == -10){ //designs less than -9
      G4double srcr,srcthet,srcz;
      srcthet=2*pi*G4UniformRand();
      srcz=(33.0)*G4UniformRand()-((33.0)/2);
      //generate uniformly in r**2/2
      srcr=(50*50/2.0)*G4UniformRand();
      srcr=sqrt(2*srcr);

      //account for the tower by adding a fixed offset here we do only the top 4 detectors of Ge!
      G4double spacing=33.0+10.16; //thickness plus 0.4in 
      G4int ndet=4;
      G4double detchoose = ndet*G4UniformRand();
      G4int thisdet=(int)detchoose+1;

      srcz += spacing*((double)thisdet-1);

      srcz += spacing/2.0;

      //set the particle position
      particleGun->SetParticlePosition(G4ThreeVector(srcr*cos(srcthet)*mm,srcr*sin(srcthet)*mm,srcz*mm).transform(*xrot));
    }
    else if(Ndet->GetDesignNo() == -11){ //designs less than -9
      G4double srcr,srcthet,srcz;
      srcthet=2*pi*G4UniformRand();
      srcz=(33.0)*G4UniformRand()-((33.0)/2);
      //generate uniformly in r**2/2
      srcr=(50*50/2.0)*G4UniformRand();
      srcr=sqrt(2*srcr);

      //account for the tower by adding a fixed offset here we do only the top 4 detectors of Ge!
      G4double spacing=33.0+10.16; //thickness plus 0.4 
      G4int ndet=4;
      G4double detchoose = ndet*G4UniformRand();
      G4int thisdet=(int)detchoose+1;

      srcz -= spacing*((double)thisdet-1);

      srcz -= spacing/2.0;

      //set the particle position
      particleGun->SetParticlePosition(G4ThreeVector(srcr*cos(srcthet)*mm,srcr*sin(srcthet)*mm,srcz*mm).transform(*xrot));
    }
    else if(Ndet->GetDesignNo() == -12){
      G4double srcr,srcthet,srcz;
      srcthet=2*pi*G4UniformRand();
      srcz=(25.4)*G4UniformRand()-((25.4)/2);
      //generate uniformly in r**2/2
      //FIXME: is this wrong? -- this is doing a diameter of 76.2 mm, so that's 3 inch
      //yep, that is the size of Soudan iZIP
      srcr=((76.2/2.0)*(76.2/2.0)/2.0)*G4UniformRand();
      srcr=sqrt(2*srcr);

      //set the particle position
      particleGun->SetParticlePosition(G4ThreeVector(srcr*cos(srcthet)*mm,srcr*sin(srcthet)*mm,srcz*mm).transform(*xrot));
    }
    else if(Ndet->GetDesignNo() == -13){
      G4double srcr,srcthet,srcz;
      srcthet=2*pi*G4UniformRand();
      srcz=(33.0)*G4UniformRand()-((33.0)/2);
      //generate uniformly in r**2/2
      srcr=(50*50/2.0)*G4UniformRand();
      srcr=sqrt(2*srcr);

      //set the particle position
      particleGun->SetParticlePosition(G4ThreeVector(srcr*cos(srcthet)*mm,srcr*sin(srcthet)*mm,srcz*mm).transform(*xrot));
    }
   
    
    else {
      G4double srcr,srcthet,srcphi;
      //generate uniform in costhet
      srcthet=-1+2*G4UniformRand();
      srcthet=acos(srcthet);
      srcphi=2*pi*G4UniformRand();
      //generate uniformly in r**3/3
      srcr=((1000.0*1000.0*1000.0)/3.0)*G4UniformRand();
      srcr=pow(3*srcr,(1/3.0));

      //set the particle position
      particleGun->SetParticlePosition(G4ThreeVector(srcr*cos(srcphi)*cos(srcthet)*mm,srcr*sin(srcphi)*cos(srcthet)*mm,srcr*sin(srcthet)*mm).transform(*xrot));
    
    }
    
    
    if(t){
      //get the energies of the cascade if the times are OK
      for(int i=0;i<ncascade;i++){
        //the decay time is actually the previous time in the chain
        double decaytime;
        if(i==0)
          decaytime=0;
        else
          decaytime=i-1;

        if(decaytime<1e12){ //1ms in fs
          angles = GenerateRandomDirection();
          particleGun->SetParticleEnergy(Eg[i]);
          particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]).transform(*xrot));
          particleGun->GeneratePrimaryVertex(anEvent);
        }
      }

      //check if the event is past the end of the file, if so, abort
      if(rootnevent>=t->GetEntries()-1){
       //hard abort
         G4cout << "about to abort" << G4endl;
         G4RunManager::GetRunManager()->AbortRun(true); 
         anEvent->SetEventAborted();
      }

      //set up for next event
      rootnevent++;
      t->GetEntry(rootnevent);
    }
    //when we want cascade but supplied no root input--generate thermal neutrons inside the material
    else{ 
      angles = GenerateRandomDirection();
      particleGun->SetParticleDefinition(G4Neutron::Definition()); 
      particleGun->SetParticleEnergy(0.0254e-6);
      particleGun->SetParticleMomentumDirection(G4ThreeVector(angles[2]*cos(angles[0]),angles[2]*sin(angles[0]),angles[1]).transform(*xrot));
      particleGun->GeneratePrimaryVertex(anEvent);

    }

  }
      
  
}
std::vector<G4double> NeutReflect_PrimaryGeneratorAction::GenerateRandomDirection()
{
  std::vector<G4double> angles;
  G4double rand1,rand2;
  rand1=G4UniformRand();
  rand2=G4UniformRand();
  angles.push_back(rand2*2*pi);            //phi
  G4double costhet = 2*rand1-1.0;
  angles.push_back(costhet);                 //costhet
  angles.push_back(sqrt(1 - costhet*costhet)); //sinthet

  return angles;

}
std::vector<G4double> NeutReflect_PrimaryGeneratorAction::GenerateGammaEnergy(G4String spectrum)
{
  std::vector<G4double> genergy;

  if(spectrum=="88Y"){
     //first generate the proper initial state distribution for 88Sr
     G4double rand1 = G4UniformRand();
     if(rand1 < 0.944){
       GenerateGammaLevel88Sr(genergy,2);
     }
     else if(rand1>=0.944 && rand1<0.999){
       GenerateGammaLevel88Sr(genergy,1);
     }
     else if(rand1>=0.999 && rand1<0.99965){
       GenerateGammaLevel88Sr(genergy,4);
     }
     else{
       GenerateGammaLevel88Sr(genergy,3);
     }

  }
  else if(spectrum=="124Sb"){
     //first generate the proper initial state distribution for 124Te 
     G4double rand1 = G4UniformRand();
     if(rand1 < 0.514){
       GenerateGammaLevel124Te(genergy,2);
     }
     else if(rand1>=0.514 && rand1<0.746){
       GenerateGammaLevel124Te(genergy,1);
     }
     else{
       GenerateGammaLevel124Te(genergy,0);
     }

  }

  return genergy;
}
void NeutReflect_PrimaryGeneratorAction::GenerateGammaLevel88Sr(std::vector<G4double> &gammas,G4int lev)
{

  if(lev==0) //base case
    return;

  G4double rand1;
  switch( lev ){

    case 1:
      gammas.push_back(1.8361);
      GenerateGammaLevel88Sr(gammas,0);
      break;

    case 2:
       rand1 = G4UniformRand();
       if(rand1 < 0.986){
         gammas.push_back(0.898);
         GenerateGammaLevel88Sr(gammas,1);
       }
       else{
         gammas.push_back(2.734);
         GenerateGammaLevel88Sr(gammas,0);
       }
       break;

    case 3:
       rand1 = G4UniformRand();
       if(rand1 < 0.75){
         gammas.push_back(1.382);
         GenerateGammaLevel88Sr(gammas,1);
       }
       else{
         gammas.push_back(3.2197);
         GenerateGammaLevel88Sr(gammas,0);
       }
       break;

    case 4:
         gammas.push_back(0.8506);
         GenerateGammaLevel88Sr(gammas,2);
	 break;
      

    default:
        break;


  }


  return;
}
void NeutReflect_PrimaryGeneratorAction::GenerateGammaLevel124Te(std::vector<G4double> &gammas,G4int lev)
{

  //FIXME this is a 2-level model, it's not great. Accounts for ~75% of decays

  if(lev==0) //base case
    return;

  G4double rand1;
  switch( lev ){

    case 1:
      gammas.push_back(0.6027);
      GenerateGammaLevel124Te(gammas,0);
      break;

    case 2:
      rand1 = G4UniformRand();
      if(rand1 < 0.928376){
        gammas.push_back(1.690971);
        GenerateGammaLevel124Te(gammas,1);
      }
      else if(rand1>=0.928376 && rand1<0.9641488){
        gammas.push_back(1.045125);
        GenerateGammaLevel124Te(gammas,0);
      }
      else if(rand1>=0.9641488 && rand1<0.99968){
        gammas.push_back(0.968195);
        GenerateGammaLevel124Te(gammas,0);
      }
      else{
        gammas.push_back(2.29402);
        GenerateGammaLevel124Te(gammas,0);
      }
      break;

    default:
        break;


  }


  return;
}
std::vector<G4double> NeutReflect_PrimaryGeneratorAction::GenerateDiskPathlength(std::vector<G4double> angles)
{
  //components of angles:
  //0: phi  1: costhet  2: sinthet

  //components of output
  //0: pathlength thru BeO
  std::vector<G4double> outpar;
 
  //intermediate vars
  G4double l,lpr,L,Lpr,thmax,thmaxpr,lstar,absd,theta;

  absd = sqrt(d*d);
  if(absd>0.0)
    thmax = atan(R/absd);
  else
    thmax = pi;
  thmaxpr = atan(R/(absd + thk)); 

  //calculate angle of real path
  if(d>0.0)
   theta = acos(angles[1]);
  else
   theta = pi - acos(angles[1]);

  //calculate some lengths
  if(theta<thmax && theta>0){
    l = tan(theta)*absd;
    lpr = std::min(R,tan(theta)*(absd+thk));
    L = l/sin(theta);
    Lpr = lpr/sin(theta);
    //lstar = sqrt((lpr-l)*(lpr-l) + thk*thk)
    lstar = (lpr-l)/sin(theta);
    outpar.push_back(lstar);
    outpar.push_back(L);
  }
  else if(theta == 0.0){
    outpar.push_back(thk);
    outpar.push_back(absd);
  }
  else{ //theta>=thmax
    outpar.push_back(0);
    outpar.push_back(0); 
  }
 
  return outpar; 
}
std::vector<G4double> NeutReflect_PrimaryGeneratorAction::GenerateShellPathlength(std::vector<G4double> angles)
{
  //components of angles:
  //0: phi  1: costhet  2: sinthet

  //components of output
  //0: pathlength thru BeO
  std::vector<G4double> outpar;

  if(d*angles[1]>0){
    outpar.push_back(thk);
    outpar.push_back(R);
  }
  else{ 
    outpar.push_back(0);
    outpar.push_back(0);
  }
 
  return outpar; 
}
G4double NeutReflect_PrimaryGeneratorAction::GenerateGammaNEnergy(G4double costhet,G4double egamma,G4double mrec, G4double Q)
{
  //get the reduced mass
  G4double mu = mrec*mn/(mrec + mn);

  //solve the quadratic equation for the 2 body --> 2 body kinematics
  G4double a = 1;
  G4double b = -egamma*mu*costhet/mrec;
  G4double c = -2*mu*(egamma + Q - (egamma*egamma)/(2*mrec));

  //G4cout << egamma << " " << mu << " " << mrec << " " << Q << " " << mn << G4endl;
  //G4cout << "Quad is solvable if: " << b*b - 4*a*c << " is greater than zero!" << G4endl;

  //quadratic formula
  G4double pn1 = (-2*b + sqrt(b*b - 4*a*c))/2;
  G4double pn2 = (-2*b - sqrt(b*b - 4*a*c))/2;

  if(pn1*pn2 > 0){
    //G4cout << "Either two roots or no solution" << G4endl;
    return 0.0; //ugh, what do I do here?
  }
  else{
    //G4cout << "One good solution!" << G4endl;
    if(pn1>0.0)
      return pn1*pn1/(2*mn);
     else
      return pn2*pn2/(2*mn);
  }

}
G4bool NeutReflect_PrimaryGeneratorAction::DidGammaNHappen(G4double pathlength,G4double energy)
{
   //initialize pathlenth
   finalpathlength = -1.0; 

   //employ threshold PR74 pg611 1949
   if(energy<1.63 || pathlength==0.0)
     return false;


   //event biasing 
   pathlength*=xnbias;
   //now that energy is above threshold calculate probability
   //this was estimated as 1.5 mb and is should be calculated as sig*t*n
   //sig is cross section, t is pathlength, n is number density
   G4double xn=1.5E-3; //in barns

   if(decayer=="88Y"){
     if(energy<=1.9)
       xn=0.748E-3;
     else if(energy>1.9 && energy<= 2.8)
       xn=0.664E-3;
     else
       xn=0.432E-3;
   }
   else if(decayer=="124Sb"){
     if(energy<=1.9)
       xn=0.73123E-3;
     else if(energy>1.9 && energy<= 2.4)
       xn=0.17575E-3;
     else
       xn=0.432E-3;
   }
   else{
     if(energy<=1.9)
       xn=0.748E-3;
     else if(energy>1.9 && energy<= 2.8)
       xn=0.664E-3;
     else
       xn=0.432E-3;
   }

   //convert to cm**2
   xn = xn*1E-24;

   //density of BeO is 3.01 g/cm3 and its molar mass is 25.01 g/mole
   //Pure Be has a molar mass of 9.012 g/mole and density of 1.85g/cm3
   G4double n;
   if(Ndet->GetBePure())
     n = 6.02E23*(1/9.012)*1.85;
   else
     n = 6.02E23*(1/25.0)*3.01;

   G4double prob = xn*(pathlength/10)*n;

   G4double rand1;
   rand1=G4UniformRand();


   if(rand1<prob){
     return true;
     //update pathlength (without bias)
     finalpathlength = pathlength/xnbias;
   }
   else
     return false;
  
}
