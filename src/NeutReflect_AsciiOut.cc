/*==================NeutReflect_AsciiOut.cc================================ 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Code support for the ascii output in the NeutReflect simulation.
              
======================================================================*/

#include "G4RunManager.hh" 
#include "G4UnitsTable.hh"

#include "NeutReflect_AsciiOut.hh"


NeutReflect_AsciiOut::NeutReflect_AsciiOut()
{

}
NeutReflect_AsciiOut::NeutReflect_AsciiOut(G4String filename, char** vnames, G4int nvar, G4int dnum)
{
  dataFileOutputName = filename;
  ///////////////  HeaderString = header;
  NumVars = nvar;
  VariableNames = vnames;
  NumDet = dnum;
  //DetNames = dnames;
  
  dataFileOutput.open(dataFileOutputName);
  dataFileOutput.precision(10);
  WriteHeaderToFile();
}
NeutReflect_AsciiOut::~NeutReflect_AsciiOut()
{
  CloseDataFile();
}
void NeutReflect_AsciiOut::DumpToFile(G4double* dataArray , G4int nrows, G4int ncols)
{
  for (G4int jj=0; jj<nrows; jj++) {
    for (G4int ii=0; ii<ncols; ii++) {
      dataFileOutput << dataArray[jj*ncols + ii] << "\t";
      //G4cout << dataArray[jj*ncols + ii] << "          " << G4endl;
    }
    dataFileOutput << G4endl;
  }  
 
}
void NeutReflect_AsciiOut::WriteHeaderToFile()
{
  for (G4int ii=0; ii<NumDet; ii++) {
    //dataFileOutput << DetNames[ii] << "\t";
  }
  dataFileOutput << G4endl;
  
  for (G4int ii=0; ii<NumVars; ii++) {
    dataFileOutput << VariableNames[ii] << "\t";
  }
  dataFileOutput << G4endl;  
}
void NeutReflect_AsciiOut::CloseDataFile()
{
  G4cout << "file name:  " << dataFileOutputName << G4endl;
  dataFileOutput.close();
}
 
// ------------------------------------------------

