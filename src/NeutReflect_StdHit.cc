/*==================NeutReflect_StdHit.cc================================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Code support for the NeutReflect_StdHit object and methods. 
              
======================================================================*/

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "NeutReflect_StdHit.hh"

G4Allocator<NeutReflect_StdHit> NeutReflect_StdHitAllocator;

NeutReflect_StdHit::NeutReflect_StdHit() 
{
  dataVector = new G4double[22];
}
NeutReflect_StdHit::NeutReflect_StdHit(const NeutReflect_StdHit& right) : G4VHit()
{
  trackID    = right.trackID;
  edep       = right.edep;
  globalTime = right.globalTime;
  pos        = right.pos;
  dataVector = right.dataVector;
}
NeutReflect_StdHit::~NeutReflect_StdHit() 
{
  delete dataVector;
}
const NeutReflect_StdHit& NeutReflect_StdHit::operator=(const NeutReflect_StdHit& right)
{
  trackID    = right.trackID;
  edep       = right.edep;
  globalTime = right.globalTime;
  pos        = right.pos;
  dataVector = right.dataVector;
  return *this;
}
int NeutReflect_StdHit::operator==(const NeutReflect_StdHit& right) const
{
  return 0;
}
void NeutReflect_StdHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
    {
      G4Circle circle(pos);
      circle.SetScreenSize(10.04);
      circle.SetFillStyle(G4Circle::filled);
      G4Colour colour(0/255.,100/255.,128/255.);
      G4VisAttributes attribs(colour);
      circle.SetVisAttributes(attribs);
      pVVisManager->Draw(circle);
    }
}
void NeutReflect_StdHit::Print()
{
  G4cout << "  trackID: " << trackID << "  Standard Sensitive" 
	 << "  energy deposit: " << G4BestUnit(edep,"Energy")
	 << "  position: " << G4BestUnit(pos,"Length") << G4endl;
}

// ------------------------------------------

