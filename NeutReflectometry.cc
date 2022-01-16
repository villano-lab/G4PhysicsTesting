/*==================NeutReflectometry.cc====================================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
 
         Anthony Villano 3/21/17:
	     Updated to be compatible with Geant4.10.01

      PURPOSE: 
              
      INPUT:      

      OUTPUT:   
              
======================================================================*/
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include "Randomize.hh"
#include <sys/time.h>
#include <time.h>
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"


#include "NeutReflect_DetectorConstruction.hh"
#include "NeutReflect_PrimaryGeneratorAction.hh"
#include "NeutReflect_RunAction.hh"
#include "NeutReflect_EventAction.hh"
#include "Shielding_ComptonsUpdate.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

G4int origevt;
G4int lineno;
// G4int thisevent;
char desiredfile[100];

//globals
long randSeed;
int main(int argc, char** argv) {

   //get the $USER env var
   char *cuser;
   cuser = getenv("USER");
   std::string user = cuser;

   G4String rep; 
   int minarg=1;
   G4String source="88Y";
   G4String nevent="10";
   G4String part="neutron";
   G4String eprimary="1.8";
   G4String strouttype = "ascii";
   G4String location="/data/vuk2/cdms/"+user+"/data/NeutReflect";
   G4String datadir="NROUT";
   G4bool nomac=false;
   G4bool incvis=false;
   G4bool UseBePure=false;
   G4bool UseDateLabel=false;
   G4int firstDesign=-1;
   G4String g4dataset="";
   G4String stringlabel="NULL";
   G4String cascaderootfile="NULL";
   G4int g4datasetno=0;
   G4double diskR=12.7;
   G4double diskt=8.0;
   G4double xnbias=1000;

   //read in arguments from command line
   for(int i=0; i<argc; i++)
   {
      rep=argv[i];
      if(rep=="-set"){
        g4datasetno=atoi(argv[i+1]);
        minarg+=2;
      }
      if(rep == "-d"){
        firstDesign=atoi(argv[i+1]);
	minarg+=2;
      }
      if(rep == "-src"){
        source = argv[i+1];
        minarg+=2;
      }
      if(rep == "-part"){
        part=argv[i+1];
	minarg+=2;
      }
      if(rep == "-E"){
        eprimary=argv[i+1];
	minarg+=2;
      }
      if(rep == "-ngen"){
        nevent=argv[i+1];
	minarg+=2;
      }
      if(rep == "-nomac"){
        nomac=true;
	minarg+=1;
      }
      if(rep == "-date"){
        UseDateLabel=true;
	minarg+=1;
      }
      if(rep == "-usestringlabel"){
	stringlabel=argv[i+1];
	minarg+=2;
      }
      if(rep == "-cascadein"){
	cascaderootfile = argv[i+1];
	minarg+=2;
      }
      if(rep == "-slac"){
        datadir="SLAC";
        location="/nfs/slac/g/cdms/u01/users/"+user+"/NeutReflect";
	minarg+=1;
      }
      if(rep == "-umn"){
        datadir="UMN";
        location="/data/chocula/"+user+"/NeutReflect";
	minarg+=1;
      }
      if(rep == "-umncascade"){
        datadir="Geant4";
        location="/data/chocula/villaa/cascadeSimData";
	minarg+=1;
      }
      if(rep=="-loc"){
        location=argv[i+1];
	minarg+=2;
      }
      if(rep=="-otype"){
        strouttype=argv[i+1];
	minarg+=2;
      }
      if(rep=="-dr"){
        diskR=atof(argv[i+1]);
	minarg+=2;
      }
      if(rep=="-dt"){
        diskt=atof(argv[i+1]);
	minarg+=2;
      }
      if(rep=="-bias"){
        xnbias=atof(argv[i+1]);
	minarg+=2;
      }
      if(rep=="-bepure"){
        UseBePure=true;
	minarg+=1;
      }
   }

   if(g4datasetno>2){
     std::ostringstream oss;
     oss << std::hex << "0x" << std::setfill('0') << std::setw(4) << g4datasetno;
     g4dataset = oss.str(); 
   }
   
   //construct location according to otype and location
   if(strouttype=="root")
     location+="/"+datadir+"root/";
   else
     location+="/"+datadir+"/";

   G4String g4design="NULL";
   if(firstDesign>=0 && firstDesign<=3){
      std::ostringstream designstream;
      designstream << "Design"<< firstDesign;
      g4design = designstream.str();
   }
   else if(firstDesign<-1){
      std::ostringstream designstream;
      designstream << "Diag"<< -firstDesign-2;
      g4design = designstream.str();

   }

   //convert radius and thickness to string
   G4String drad,dthick,bias;
   std::ostringstream str_drad,str_dthick,str_bias;
   str_drad << diskR;
   str_dthick << diskt;
   str_bias << xnbias;
   drad = str_drad.str();
   dthick = str_dthick.str();
   bias = str_bias.str();
   G4cout << drad << G4endl;
   G4cout << dthick << G4endl;
   G4cout << bias << G4endl;

   if(drad.find(".",0)!=std::string::npos)
     drad.replace(drad.find(".",0),1,"-");
   if(dthick.find(".",0)!=std::string::npos)
     dthick.replace(dthick.find(".",0),1,"-");
   if(bias.find(".",0)!=std::string::npos)
     bias.replace(bias.find(".",0),1,"-");

#ifdef G4VIS_USE
  // Visualisation, if you choose to have it!
  incvis=true;
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  //Run Manager
  G4RunManager * runManager = new G4RunManager;

  // User Initializaton classes (mandatory)
  runManager->SetUserInitialization(new NeutReflect_DetectorConstruction(firstDesign,diskR,diskt,UseBePure));
  //runManager->SetUserInitialization(new NeutReflect_PhysicsList);
  //runManager->SetUserInitialization(new QGSP_BERT_HP);
  //runManager->SetUserInitialization(new CDMS_SuperSim_PhysicsList);
  runManager->SetUserInitialization(new Shielding_ComptonsUpdate);

  // set manditory user action class
  // UserAction Classes
  NeutReflect_PrimaryGeneratorAction* myPrimaryEventGenerator;
  if(cascaderootfile=="NULL")
    myPrimaryEventGenerator=new NeutReflect_PrimaryGeneratorAction(xnbias,source);
  else
    myPrimaryEventGenerator=new NeutReflect_PrimaryGeneratorAction(xnbias,"cascade",cascaderootfile);
  runManager->SetUserAction(myPrimaryEventGenerator);
  NeutReflect_RunAction* pRunAction = new NeutReflect_RunAction(myPrimaryEventGenerator);

  runManager->SetUserAction(pRunAction);
  runManager->SetUserAction(new NeutReflect_EventAction(pRunAction));

  // Initialize G4 kernel
  runManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager * UI = G4UImanager::GetUIpointer();

  //create a session just in case
  G4UIsession *session = 0;
  session = new G4UIterminal();

    if(argc==minarg && !nomac)
    // Define (G)UI terminal for interactive mode
    {
      if(!incvis){
        G4cerr << "GeantFlukaCompare: ERROR! code was compiled without visualization/graphics support"
	  << G4endl;
        G4cerr << "GeantFlukaCompare: default.mac includes graphics so cannot be run in this mode"
	  << G4endl;
	exit(1);
      }
      G4String command;
      G4String modifier="";
      if(diskR!=12.7 || diskt!=8.0){
        modifier="BeR_"+drad+"_Bet_"+dthick;
        if(UseBePure)
          modifier += "BePure";
        //I guess 1000 has become my standard bias
        if(xnbias!=1000){
          modifier += "bias_"+bias;
        }
      }
      G4String fileprefix;
      if(UseDateLabel || stringlabel!="NULL"){

	//turn off random numbers
        command = "/run/NeutReflect/UseRandomTag false";
        UI->ApplyCommand(command);

	//either get a time or a stringlabel bassed on the inputs
        if(stringlabel!="NULL") //let the string label take presedence
          fileprefix = location+"NeutReflect"+g4design+g4dataset+"_Source"+source+"_"+part+"_"+stringlabel;
	else{
          //FIXME should do date but I'm too lazy right now
          fileprefix = location+"NeutReflect"+g4design+g4dataset+"_Source"+source+"_"+part;
	}
      }
      else
        fileprefix = location+"NeutReflect"+g4design+g4dataset+"_Source"+source+"_"+part;
      fileprefix = fileprefix+modifier;
      command = "/run/NeutReflect/OFPrefix "+fileprefix;
      UI->ApplyCommand(command);

      UI->ApplyCommand("/control/execute default.mac");
      session->SessionStart();
      delete session;
    }
    else if(argc==minarg && nomac)
    // Define (G)UI terminal for interactive mode
    {
      G4String g4part="gamma";
      if(part=="mu")
        g4part="mu-";
      else if(part=="amu")
        g4part="mu+";
      else
        g4part=part;


      G4String command;
      G4String modifier="";
      if(diskR!=12.7 || diskt!=8.0){
        modifier="BeR_"+drad+"_Bet_"+dthick;
        if(UseBePure)
          modifier += "BePure";
        //I guess 1000 has become my standard bias
        if(xnbias!=1000){
          modifier += "bias_"+bias;
        }
      }
      G4String fileprefix;
      if(UseDateLabel || stringlabel!="NULL"){

	//turn off random numbers
        command = "/run/NeutReflect/UseRandomTag false";
        UI->ApplyCommand(command);

	//either get a time or a stringlabel bassed on the inputs
        if(stringlabel!="NULL") //let the string label take presedence
          fileprefix = location+"NeutReflect"+g4design+g4dataset+"_Source"+source+"_"+part+"_"+stringlabel;
	else{
          //FIXME should do date but I'm too lazy right now
          fileprefix = location+"NeutReflect"+g4design+g4dataset+"_Source"+source+"_"+part;
	}
      }
      else
        fileprefix = location+"NeutReflect"+g4design+g4dataset+"_Source"+source+"_"+part;
      
      fileprefix = fileprefix+modifier;
      command = "/run/NeutReflect/OFPrefix "+fileprefix;
      UI->ApplyCommand(command);
      //command = "/run/NeutReflect/OutType "+strouttype;
      //UI->ApplyCommand(command);
      command = "/run/beamOn "+nevent;
      UI->ApplyCommand(command);
    }
  else
    // Batch mode
    {
      if(!incvis){
        G4cout << "GeantFlukaCompare: WARNING! code was compiled without visualization/graphics support"
	  << G4endl;
      }
      G4String command = "/control/execute ";
      G4String macroFileName = argv[minarg];
      G4cout << " Macro " << macroFileName << " is being run " << G4endl;
      UI->ApplyCommand(command+macroFileName);
    
      session->SessionStart();
      delete session;
    }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
