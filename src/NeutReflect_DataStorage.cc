/*==================NeutReflect_DataStorage.cc============================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12

      UPDATES:      
       

      PURPOSE:  Code support for the data storage class and routines
                to handle the organization output data into various
		forms which default to ascii supported by NeutReflect_AsciiOut.
              
======================================================================*/

#include "NeutReflect_DataStorage.hh"

#ifdef ASCIIOUT
#include "NeutReflect_AsciiOut.hh"
#endif

#include "NeutReflect_DetectorConstruction.hh"

NeutReflect_DataStorage::NeutReflect_DataStorage()
{
  printf(" NeutReflect_DATA storage -- EMPTY constructor\n");
}
NeutReflect_DataStorage::NeutReflect_DataStorage(G4String filename, G4int run, G4int rseed)
{  

  NumDets = 1;

  //  storeEvtPrimaries = 0; //// This should come from RunAction (which is set via RunMessenger)

  n_data    = N_DATA;
  n_length  = N_LENGTH;
  n_entries = 0;
  n_files   = 0;

  runID = run;
  randSeed = rseed;
  outfilename = filename;


  data = new G4double [N_DATA];
  dataArray = new G4double[N_LENGTH*N_DATA];
 
  
  VariableNames = new char*[N_DATA];
  char myNames[N_DATA][20];
  sprintf(myNames[0], "EV"); sprintf(myNames[1], "DT"); sprintf(myNames[2], "TS");
  sprintf(myNames[3], "P"); sprintf(myNames[4], "Type"); sprintf(myNames[5], "E1");
  sprintf(myNames[6], "D3");
  sprintf(myNames[7], "X1"); sprintf(myNames[8], "Y1"); sprintf(myNames[9], "Z1");
  sprintf(myNames[10], "X3"); sprintf(myNames[11], "Y3"); sprintf(myNames[12], "Z3");
  sprintf(myNames[13], "PX"); sprintf(myNames[14], "PY"); sprintf(myNames[15], "PZ");
  sprintf(myNames[16], "time");
  sprintf(myNames[17], "origevt");

  for (int kk=0; kk<N_DATA; kk++)
    VariableNames[kk] = myNames[kk];


#ifdef ASCIIOUT

  // Open the data file for writing 
  G4cout << "DataStorage sees filename as: " << filename << " upon construction" << G4endl;
  outfilename = filename + G4String("_") + NumtoStr(runID,3) + G4String("_"); 
  Out = new NeutReflect_AsciiOut(outfilename + NumtoStr(n_files,3) + ".txt", VariableNames, n_data, NumDets);

#endif

}
NeutReflect_DataStorage::~NeutReflect_DataStorage()
{
#ifdef ASCIIOUT
  Out->DumpToFile(dataArray , n_entries, n_data);
  delete Out;
#endif

  //delete data;  
  delete dataArray;
  G4cout << "Deleted dataArray." << G4endl;
}
void NeutReflect_DataStorage::setData(G4double* input)
{
  data = input;

}
G4double* NeutReflect_DataStorage::getData()
{
  return data;
}
void NeutReflect_DataStorage::addData()
{
  //  for (G4int ii=0; ii<N_DATA; ii++) {
  //    dataArray[n_entries][ii] = data[ii];
  //  }
  //G4cout << "n_data " << n_data << "  and n_entries=  " << n_entries << G4endl;
 
  if (n_entries>n_length)
    writeArray();

  for (G4int ii=0; ii<n_data; ii++) {
    dataArray[n_entries*n_data + ii] = data[ii];
  }
  n_entries++;  



}
void NeutReflect_DataStorage::resetArray()
{

  for (G4int jj=0; jj<n_entries; jj++) {
    for (G4int ii=0; ii<n_data; ii++) {
      dataArray[jj*n_data + ii] = -99999;
      }
  }
  n_entries = 0;  
}
void NeutReflect_DataStorage::writeArray()
{

  char myNames[N_DATA][20];
  sprintf(myNames[0], "EV"); sprintf(myNames[1], "DT"); sprintf(myNames[2], "TS");
  sprintf(myNames[3], "P"); sprintf(myNames[4], "Type"); sprintf(myNames[5], "E1");
  sprintf(myNames[6], "D3");
  sprintf(myNames[7], "X1"); sprintf(myNames[8], "Y1"); sprintf(myNames[9], "Z1");
  sprintf(myNames[10], "X3"); sprintf(myNames[11], "Y3"); sprintf(myNames[12], "Z3");
  sprintf(myNames[13], "PX"); sprintf(myNames[14], "PY"); sprintf(myNames[15], "PZ");
  sprintf(myNames[16], "time");
  sprintf(myNames[17], "origevt");

  for (int kk=0; kk<N_DATA; kk++)
    VariableNames[kk] = myNames[kk];



#ifdef ASCIIOUT
  Out->DumpToFile(dataArray , n_entries, n_data);
  delete Out;
  n_files++;

  Out = new NeutReflect_AsciiOut(outfilename + NumtoStr(n_files,3) + G4String(".txt"), VariableNames, n_data, NumDets);
#endif

  resetArray();
}
void NeutReflect_DataStorage::printArray()
{

  for (G4int jj=0; jj<n_entries; jj++) {
    G4cout << "******** DataStorage : -----" << G4endl;
    for (G4int ii=0; ii<n_data; ii++) {
      G4cout << "data["<< jj << "," << ii <<  "] = " << dataArray[jj*n_data + ii] << " - ";
      }
    G4cout << G4endl;
  }
  G4cout << G4endl;

}
G4String NeutReflect_DataStorage::NumtoStr(G4int number, G4int length)
{
  //G4String name;
  std::ostringstream oss;

  oss << std::setfill('0');
  oss << std::setw(length) << number;

  return oss.str();
}
G4bool NeutReflect_DataStorage::overflowArray(G4int hits)
{
  if (n_entries+hits > n_length) return TRUE;
  else return FALSE;
}

