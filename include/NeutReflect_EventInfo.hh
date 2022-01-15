/*==================NeutReflect_EventInfo.hh========================= 
   
      PROGRAMMER:  Anthony Villano  04/10/15

      UPDATES:      
       

      PURPOSE: Class for a NeutReflect to store some add'l event info
              
======================================================================*/

#ifndef NeutReflect_EventInfo_h
#define NeutReflect_EventInfo_h 1

// ------------------------------------------------

#include "globals.hh"


// ------------------------------------------------

class NeutReflect_EventInfo : public G4VUserEventInformation
{
public :
  NeutReflect_EventInfo();
  ~NeutReflect_EventInfo();

  //pure virtuals
  void Print() const;

  G4double GetBeLength(){return BeLength;}
  void SetBeLength(G4double length){BeLength = length;}

public :
  
  
private :

  G4double BeLength;

};

#endif
