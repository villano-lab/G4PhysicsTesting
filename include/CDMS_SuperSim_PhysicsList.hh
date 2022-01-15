#ifndef CDMS_SuperSim_PhysicsList_hh
#define CDMS_SuperSim_PhysicsList_hh 1
// $Id: CDMS_SuperSim_PhysicsList.hh,v 1.1 2012/10/20 00:47:32 villaa Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSPhysicsList.hh                                   //
//  Description: Physics processes relevant to CDMS backgrounds       //
//                                                                    //
//  Author:      Michael Kelsey                                       //
//  Date:        14 February 2011                                     //
//                                                                    //
//  20110322  M. Kelsey -- Add optical-photon physics.                //
//  20110427  M. Kelsey -- Add verbose as constructor argument.       //
//  20110525  M. Kelsey -- Add lists' instantiation function.         //
//  20111129  M. Kelsey -- Drop data members; inherit from Shielding  //
////////////////////////////////////////////////////////////////////////

#include "Shielding.hh"
#include "globals.hh"

class CDMS_SuperSim_PhysicsMessenger;

// FIXME:  Can't use typedef name for inheritance?
class CDMS_SuperSim_PhysicsList: public Shielding {
public:
  CDMS_SuperSim_PhysicsList(G4int verbose=1);	// Physics lists use "1" for quiet
  virtual ~CDMS_SuperSim_PhysicsList();

private:
  CDMS_SuperSim_PhysicsMessenger* messenger;
};

#endif	/* CDMS_SuperSim_PhysicsList_hh */
      
