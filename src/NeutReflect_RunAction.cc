/*==================NeutReflect_RunAction.cc============================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE: Code supporting the RunAction class including specific
               routines for setting file name prefix and beinginning/end
	       of run actions.
              
======================================================================*/

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

#include "Randomize.hh"

#include "NeutReflect_RunAction.hh"
#include "NeutReflect_RunActionMessenger.hh"
#include "NeutReflect_DetectorConstruction.hh"
#include "NeutReflect_DataStorage.hh"
#include "NeutReflect_PrimaryGeneratorAction.hh"
#include <sys/time.h>
#include <time.h>
char filename[200];
// extern G4int thisevent;


NeutReflect_RunAction::NeutReflect_RunAction(NeutReflect_PrimaryGeneratorAction* generator)
{
  
  
  // automatic (time-based) random seeds and filenames for each run
  struct timeval mytime;
  gettimeofday(&mytime, NULL);
  randSeed=mytime.tv_sec-mytime.tv_usec;
  CLHEP::HepRandom::setTheSeed(randSeed);
  CLHEP::HepRandom::showEngineStatus();

  // Set a default file prefix
  sprintf(filename,"NeutReflect_%ld",randSeed);
  DataFileNamePrefix = G4String(filename);
  UseRandomTag=true;
  OutputDataToFile = true;
  DrawEventCmd = true;
  SaltPillOutCmd = false;
  thisgenerator=generator;

  runMessenger = new NeutReflect_RunActionMessenger(this);  
}
NeutReflect_RunAction::~NeutReflect_RunAction()
{

}
void NeutReflect_RunAction::BeginOfRunAction(const G4Run* aRun)
{
  // thisevent=0;

  runN = aRun->GetRunID();
  if ( (runN % 1000 == 0) || (runN<100) ) 
    G4cout << "### Run : " << runN << G4endl;

  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer();
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    }

  // ------ Determine data output file name -----

  // There MUST be a better way of doing this !!!
  char fnum[5];
  sprintf(fnum,"%03d",runN);

  if(OutputDataToFile) {

    // Create the DataOut instance -- with the current data file names
    char file[200];
    sprintf(file,"%ld",randSeed);
    G4cout << "Prior to constructing DataStorage RunAction sees DataFileNamePrefix as: " << DataFileNamePrefix << G4endl;
    if(false) //FIXME just hack this out for now
      dataOut = new NeutReflect_DataStorage(DataFileNamePrefix+G4String(file),runN,1); 
    else
      dataOut = new NeutReflect_DataStorage(DataFileNamePrefix,runN,1); 

  }

}
void NeutReflect_RunAction::EndOfRunAction(const G4Run*)
{

  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");

  //Kill the data out class instance
  if(OutputDataToFile) {
    delete dataOut;
  }

}
