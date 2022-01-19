This is a command that works on my mac:

* create a directory in the build directory called NROUT

* run: ./NeutReflectometry -src PuBe -d -1 -set 0 -ngen 100000 -otype txt -usestringlabel test -nomac -loc .

* the  -d option selects the detector design, essentially the geometry. Here is a list of choices:

```
	
	//------------------------------------------------ 
	// Construct other detector components 
	//------------------------------------------------ 
        //Parameters
        G4double dist=400.0;
        G4double tl=6.0,tu=100.0,dr=30.0,g1=10.0,g2=10.0,hr=200.0,rr=30.0;
        G4double dpipe=5.0,tsteelplate=20.0,tsrcholder=6.35,lBeO=400.0,ldet=186.2,ddet=152.4;
        //G4double dpipe=5.0,tsteelplate=6.35,tsrcholder=6.35,lBeO=400.0,ldet=186.2,ddet=152.4;
        //start by constructing Ge detector
        if(DesignNo>-2){ //use below -2 for diagnostic geometries
          ConstructGeDet(physicalWorld,dist);
          ConstructSimpleVessel(physicalWorld,40.0,76.2+dr+g2+rr);
        }

        if(DesignNo==0){
          //construct Design0
          ConstructDesign0(physicalWorld,tl);
        }
        else if(DesignNo==1){
          //construct Design1
          ConstructDesign1(physicalWorld,tl,tu,dr,g1,g2,hr,rr);
        }
        else if(DesignNo==2){
          //construct Design2
          ConstructDesign2(physicalWorld,tl,dist);
        }
        else if(DesignNo==3){
          //construct Design3
          ConstructDesign3(physicalWorld,tl,dpipe,tsteelplate,tsrcholder,lBeO,ldet,ddet);
        }
        else if(DesignNo==-1){
          //construct NULL
          ConstructNULL(physicalWorld);
        }
        else if(DesignNo==-2){
          //construct Diag0 
          ConstructDiag0(physicalWorld);
        }
        else if(DesignNo==-3){
          //construct Diag1 
          ConstructDiag1(physicalWorld);
        }
        else if(DesignNo==-4){
          //construct Diag2 
          ConstructDiag2(physicalWorld);
        }
        else if(DesignNo==-5){
          //construct Diag3
	  //same as Diag2 but with no lead stack
          ConstructDiag3(physicalWorld);
        }
        else if(DesignNo==-6){
          //construct Diag4
	  //same as Diag3 but with nothing but Be disk
          ConstructDiag4(physicalWorld);
        }
        else if(DesignNo==-7){
          //construct Diag5
	  //just a Ge detector centered at zero
          ConstructGeDet(physicalWorld,0.0);
        }
        else if(DesignNo==-8){
          //construct Diag6
	  //just a Si detector centered at zero
          ConstructSiDet(physicalWorld,0.0);
        }
        else if(DesignNo==-9){
          //construct Diag7
	  //just a Ar detector centered at zero
          ConstructArDet(physicalWorld,0.0);
        }
        else if(DesignNo==-10){
          //construct Diag8
	  //just a Ge/Si tower of SNOLAB dets 
	  //for cascade generation in the Ge detectors ONLY
	  G4double spacing = 33.0+10.16; //use 0.4in spacing
	  //spacing negative switches order of towers
	  spacing *= -1;
          ConstructSNOLABDet(physicalWorld,"Ge1","Ge",(spacing/2.0)+3*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge2","Ge",(spacing/2.0)+2*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge3","Ge",(spacing/2.0)+1*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge4","Ge",(spacing/2.0)+0*spacing);
          ConstructSNOLABDet(physicalWorld,"Si5","Si",-(spacing/2.0)-0*spacing);
          ConstructSNOLABDet(physicalWorld,"Si6","Si",-(spacing/2.0)-1*spacing);
          ConstructSNOLABDet(physicalWorld,"Si7","Si",-(spacing/2.0)-2*spacing);
          ConstructSNOLABDet(physicalWorld,"Si8","Si",-(spacing/2.0)-3*spacing);
        }
        else if(DesignNo==-11){
          //construct Diag9
	  //just a Ge/Si tower of SNOLAB dets 
	  //for cascade generation in the Si detectors ONLY
	  G4double spacing = 33.0+10.16; //use 0.4in spacing
	  //spacing negative switches order of towers
	  spacing *= -1;
          ConstructSNOLABDet(physicalWorld,"Ge1","Ge",(spacing/2.0)+3*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge2","Ge",(spacing/2.0)+2*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge3","Ge",(spacing/2.0)+1*spacing);
          ConstructSNOLABDet(physicalWorld,"Ge4","Ge",(spacing/2.0)+0*spacing);
          ConstructSNOLABDet(physicalWorld,"Si5","Si",-(spacing/2.0)-0*spacing);
          ConstructSNOLABDet(physicalWorld,"Si6","Si",-(spacing/2.0)-1*spacing);
          ConstructSNOLABDet(physicalWorld,"Si7","Si",-(spacing/2.0)-2*spacing);
          ConstructSNOLABDet(physicalWorld,"Si8","Si",-(spacing/2.0)-3*spacing);
        }
        else if(DesignNo==-12){
          //construct Diag10
	  //just a Ge detector centered at zero
          ConstructGeDetSoudan(physicalWorld,0.0);
        }
        else if(DesignNo==-13){
	  //just a Si detector centered at zero
          ConstructR68SiDet(physicalWorld,0.0);
        }

```
