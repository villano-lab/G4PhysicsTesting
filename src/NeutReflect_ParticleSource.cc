/*==================NeutReflect_ParticleSource.cc========================== 
   
      PROGRAMMER:  Anthony Villano 10/18/12 

      UPDATES:      
       

      PURPOSE: Code support for raA particle source similar to the CDMS_ParticleSource.cc
               used in the /cdmsim/Gopher/ simulation.  The specialty of
	       this type of source is generating contamination by generating
	       all the gamma sources which come from gamma decays.  The gammas
	       are hard-coded and so are the event generators which generate
	       throughout the volumes.  It is necessary to change this file
	       directly if either the source changes or the geometry changes.
	       This is a major shortfall of this module, and it should therefore
	       be made to use the geometry and possibly a source input file to
	       look up sources in a database and also dynamically generate in 
	       geometry which actually exists in a simulation (tagging on a
	       material would be phenomenal). 
              
======================================================================*/
