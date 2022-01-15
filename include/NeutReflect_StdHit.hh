/*==================NeutReflect_StdHit.hh================================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Class for specifying a sensitive detector hits collections
               in the NeutReflect simulation. 
              
======================================================================*/

#ifndef NeutReflect_StdHit_h
#define NeutReflect_StdHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

// ------------------------------------------

class NeutReflect_StdHit : public G4VHit
{
public:
  NeutReflect_StdHit();
  NeutReflect_StdHit(const NeutReflect_StdHit&);
	
  ~NeutReflect_StdHit();
	
  const NeutReflect_StdHit& operator=(const NeutReflect_StdHit&);
  int operator==(const NeutReflect_StdHit&) const;
	
  inline void* operator new(size_t);
  inline void  operator delete(void*);
	
  void Draw();
  void Print();
	
public:
  void SetTrackID  (G4int track)       {trackID = track;};
  void SetDetNb    (G4int det)         {detNb = det;};
  void SetGlobalTime(G4double gTime) {globalTime = gTime;};
  void SetEdep     (G4double de)       {edep = de;};
  void SetPos      (G4ThreeVector xyz) {pos = xyz;};
  void SetPName    (G4String pName)    {particleName = pName;};
	
  G4int     GetTrackID(){return trackID;};
  G4int       GetDetNb(){return detNb;};
  G4double     GetEdep(){return edep;};
  G4double GetGlobalTime() {return globalTime;};
  G4ThreeVector GetPos(){return pos;};
  G4String GetParticle(){return particleName;};
  G4double*    GetData(){return dataVector;};

private:
  G4int         trackID;
  G4int         detNb;
  G4double      edep;
  G4double      globalTime;
  G4ThreeVector pos;
  G4String      particleName;
  G4double*     dataVector;
};

typedef G4THitsCollection<NeutReflect_StdHit> NeutReflect_StdHitsCollection;

extern G4Allocator<NeutReflect_StdHit> NeutReflect_StdHitAllocator;

inline void* NeutReflect_StdHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) NeutReflect_StdHitAllocator.MallocSingle();
  return aHit;
}
inline void NeutReflect_StdHit::operator delete(void *aHit)
{
  NeutReflect_StdHitAllocator.FreeSingle((NeutReflect_StdHit*) aHit);
}

#endif

