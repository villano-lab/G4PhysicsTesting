/*==================NeutReflect_DataStorage.hh============================ 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE: Class to specify the format for data storage and handle
               calls to the NeutReflect_AsciiOut class, which writes the
	       data.
              
======================================================================*/

#ifndef NeutReflect_DataStorage_H
#define NeutReflect_DataStorage_H 1

#define ASCIIOUT 1


#include "globals.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"

#include <list>

#define N_DATA 18
#define N_LENGTH 200000

#ifdef ASCIIOUT
class NeutReflect_AsciiOut;
#endif

class NeutReflect_DataStorage
{
public:
  NeutReflect_DataStorage();
  NeutReflect_DataStorage(G4String, G4int, G4int);
	
  ~NeutReflect_DataStorage();
	
private:
  G4int storeEvtPrimaries,lastGE;
  G4int n_data,n_length; 
  G4int n_adds,n_entries,n_files;
  
  G4int runID, randSeed;
  G4String outfilename;

  G4double* data;
  G4double* dataArray;

  char** DetNames;
  char** VariableNames;
  G4int NumDets;

#ifdef ASCIIOUT
  NeutReflect_AsciiOut* Out;
#endif

public:		
  // member functions
  void setData(double*);
  G4double* getData();
  void addData();
  void printArray();
  void writeArray();
  G4bool overflowArray(G4int);

private:
  void resetArray();
  G4String NumtoStr(G4int, G4int);

public:

};

#endif
