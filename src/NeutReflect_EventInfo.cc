/*==================NeutReflect_EventInfo.cc============================ 
   
      PROGRAMMER:  Anthony Villano 04/10/15 

      UPDATES:      
       

      PURPOSE: Code support for the ascii output in the NeutReflect simulation.
              
======================================================================*/

#include "G4RunManager.hh" 
#include "G4UnitsTable.hh"

#include "NeutReflect_EventInfo.hh"


NeutReflect_EventInfo::NeutReflect_EventInfo():
BeLength(-1.0)
{

}
NeutReflect_EventInfo::~NeutReflect_EventInfo()
{
  
}
void NeutReflect_EventInfo::Print() const
{
  G4cout << "The gamma path length through the Be: " << BeLength << G4endl;

  return;
}
// ------------------------------------------------

