/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.0 $// $Date: 2014-10-20 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/MultipleShearSpring_forThermalConductivityAnalysis/MultipleShearSpring_forThermalConductivityAnalysis.h,v $

// Written: Ippei Tsunezawa
// Created: October 2014
//
// Multiple Shear Spring (MSS) model for Thermal Conductivity Analysis
//
// Description: This file contains the implementation of the MultipleShearSpring_forThermalConductivityAnalysis class.
//              This file contains the function to parse the TCL input
// element MultipleShearSpring_forThermalConductivityAnalysis eleTag? iNode? jNode? nSpring? -mat matTag? <-lim limDisp?> <-dtstep dt?> <-mesh tate? yoko? hyoko? ctate? ftate? jtate? itate? norl? noss?> <-height hc? hf? hj? hrl? hss? hi?> <-temptarget targetorauto? targettate? targetyoko?> <-loop l?> <-totalnumber n?> <-physicalproperty capap? capaf? capac? caparl? condp? condf? condc? condrl> <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>

// This element works in 3-dim 6-dof model.
// The multiple shear springs are distributed in local-yz plane.
// Stiffness of local-x,rx,ry,rz remain zero.

#include <MultipleShearSpring_forThermalConductivityAnalysis.h>
#include <TclModelBuilder.h>

#include <Domain.h>
#include <ID.h>
#include <Vector.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <KikuchiAikenLRB_TCA.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>


extern void printCommand(int argc, TCL_Char **argv);


int TclModelBuilder_addMultipleShearSpring_forThermalConductivityAnalysis(ClientData clientData,
    Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
    TclModelBuilder *theTclBuilder)
{

  // ensure the destructor has not been called
  if (theTclBuilder == 0)  {
    opserr << "WARNING builder has been destroyed - MultipleShearSpring_forThermalConductivityAnalysis\n";
    return TCL_ERROR;
  }

  // 3-dim, 6-dof
  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING MultipleShearSpring_forThermalConductivityAnalysis command only works when ndm is 3 and ndf is 6" << endln;
    return TCL_ERROR;
  }



  //arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;
  int nSpring;
  int matTag;

  //material
  KikuchiAikenLRB_TCA *material;
  int recvMat = 0;

  //arguments (optional)
  double limDisp = 0.0;
  Vector oriX(0);
  Vector oriYp(3); oriYp(0) = 0.0; oriYp(1) = 1.0; oriYp(2) = 0.0;
  double mass = 0.0;


  //デフォルト値の設定
  int    tate    = 36;
  int    yoko    = 25;
  int    hyoko   = 4;
  int    ctate   = 1;
  int    ftate   = 1;
  int    jtate   = 1;
  int    itate   = 0;
  int    norl    = 33;
  int    noss    = 32;
  double hc      = 0.050;
  double hf      = 0.030;
  double hj      = 0.035;
  double hrl     = 0.003;
  double hss     = 0.0022;
  double hi      = 0.0;
  double dtstep  = 0.01;
   //double dtstep = ops_Dt;
  int    targettate   = 36;
  int    targetyoko   = 25;
  int    loop         = 5;
  int    totalnumber  = 1;
  int    target_type  = 0;
  double capap        = 1472900.0;
  double capaf        = 3782400.0;
  double capac        = 2409600.0;
  double caparl       = 1753050.0;
  double condp        = 35.2;
  double condf        = 45.0;
  double condc        = 2.1;
  double condrl       = 1.42;


  Element *theElement = 0;


  //error flag
  bool ifNoError = true;


  if (argc < 8)  { //element MultipleShearSpring_forThermalConductivityAnalysis eleTag? iNode? jNode? nSpring? -mat matTag?

    opserr << "WARNING insufficient arguments\n";
    ifNoError = false;

  } else {

    //argv[2~5]
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK)  {
      opserr << "WARNING invalid MultipleShearSpring_forThermalConductivityAnalysis eleTag\n";
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK)  {
      opserr << "WARNING invalid iNode\n";
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK)  {
      opserr << "WARNING invalid jNode\n";
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[5], &nSpring) != TCL_OK || nSpring <= 0)  {
      opserr << "WARNING invalid nSpring\n";
      ifNoError = false;
    }


    //argv[6~]
    for (int i=6; i<=(argc-1); i++) {
      
      double value;
      
      if (strcmp(argv[i],"-mat")==0 && (i+1)<=(argc-1)) { // -mat matTag?
	
	if (Tcl_GetInt(interp,argv[i+1], &matTag) != TCL_OK) {
	  opserr << "WARNING invalid matTag\n";
	  ifNoError = false;
	}


	material = (KikuchiAikenLRB_TCA *)OPS_getUniaxialMaterial(matTag);
	if (material == 0)  {
	  opserr << "WARNING material model not found\n";
	  opserr << "uniaxialMaterial: " << matTag << endln;
	  opserr << "MultipleShearSpring_forThermalConductivityAnalysis element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	//opserr << "org material " << material->getClassType() << "\n";
	recvMat++ ;
	i += 1;

      } else if (strcmp(argv[i],"-orient")==0 && (i+6)<=(argc-1) && Tcl_GetDouble(interp,argv[i+4], &value) == TCL_OK) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

	oriX.resize(3);

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    opserr << "WARNING invalid -orient value\n";
	    ifNoError = false;
	  } else {
	    oriX(j-1) = value;
	  }
	}

	i += 3;

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    opserr << "WARNING invalid -orient value\n";
	    ifNoError = false;
	  } else {
	    oriYp(j-1) = value;
	  }
	}

	i += 3;


      } else if (strcmp(argv[i],"-orient")==0 && (i+3)<=(argc-1)) { // <-orient yp1? yp2? yp3?> の読み込み

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    opserr << "WARNING invalid -orient value\n";
	    ifNoError = false;
	  } else {
	    oriYp(j-1) = value;
	  }
	}


	i += 3;

      } else if (strcmp(argv[i],"-mass")==0 && (i+1)<=(argc-1)) { // <-mass m?> の読み込み
	
	if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK || mass <= 0)  {
	  opserr << "WARNING invalid mass\n";
	  ifNoError = false;
	}

	i += 1;

      } else if (strcmp(argv[i],"-dtstep")==0 && (i+1)<=(argc-1)) { // <-dtstep dt?> の読み込み
	
	if (Tcl_GetDouble(interp, argv[i+1], &dtstep) != TCL_OK || dtstep <= 0)  {
	  opserr << "WARNING invalid dtstep\n";
	  ifNoError = false;
	}

	i += 1;


      } else if (strcmp(argv[i],"-mesh")==0 && (i+9)<=(argc-1)) { // <-mesh tate? yoko? hyoko? ctate? ftate? jtate? itate? norl? noss?> の読み込み

	if (Tcl_GetInt(interp,argv[i+1], &tate) != TCL_OK || tate <= 0 ) {
	  opserr << "WARNING invalid tate" << endln;
	  ifNoError = false;
	}

	if (Tcl_GetInt(interp,argv[i+2], &yoko) != TCL_OK || yoko <= 0) {
	  opserr << "WARNING invalid yoko" << endln;
	  ifNoError = false;
	}

	if (Tcl_GetInt(interp, argv[i+3], &hyoko) != TCL_OK )  {
	  opserr << "WARNING invalid hyoko\n";
	  ifNoError = false;
	}

	if (Tcl_GetInt(interp,argv[i+4], &ctate) != TCL_OK || ctate <= 0) {
	  opserr << "WARNING invalid ctate" << endln;
	  ifNoError = false;
	}

	if (Tcl_GetInt(interp,argv[i+5], &ftate) != TCL_OK || ftate <= 0) {
	  opserr << "WARNING invalid ftate" << endln;
	  ifNoError = false;
	}

	if (Tcl_GetInt(interp, argv[i+6], &jtate) != TCL_OK )  {
	  opserr << "WARNING invalid jtate\n";
	  ifNoError = false;
	}

	if (Tcl_GetInt(interp, argv[i+7], &itate) != TCL_OK )  {
	  opserr << "WARNING invalid itate\n";
	  ifNoError = false;
	}

	if (Tcl_GetInt(interp, argv[i+8], &norl) != TCL_OK )  {
	  opserr << "WARNING invalid norl\n";
	  ifNoError = false;
	}

	if (Tcl_GetInt(interp, argv[i+9], &noss) != TCL_OK )  {
	  opserr << "WARNING invalid noss\n";
	  ifNoError = false;
	}

	i += 9;


      } else if (strcmp(argv[i],"-temptarget")==0 && (i+3)<=(argc-1) && Tcl_GetInt(interp,argv[i+2], &targettate) == TCL_OK) { //<-temptarget targetorauto? targettate? targetyoko?> の読み込み
 
	if (( strcmp(argv[i+1],"target") == 0)) {
	  target_type=0;//固定のメッシュで温度を採用する際はtarget_type=0
	}else{
	  opserr << "WARNING invalid temptarget" << endln;
	  ifNoError = false;
	}
      
	if (Tcl_GetInt(interp, argv[i+2], &targettate) != TCL_OK || targettate <= 0)  {
      	  opserr << "WARNING invalid targettate\n";
      	  ifNoError = false;
      	}

      	if (Tcl_GetInt(interp, argv[i+3], &targetyoko) != TCL_OK || targetyoko <= 0)  {
      	  opserr << "WARNING invalid targetyoko\n";
      	  ifNoError = false;
      	}



      	i += 3;

      } else if (strcmp(argv[i],"-temptarget")==0 && (i+1)<=(argc-1) ) { //鉛プラグ温度の最大値を採用温度(temptarget)とする際の読み込み
       
	if (( strcmp(argv[i+1],"auto") == 0)) {
	  target_type=1;//最大温度を引き出す際はtarget_type=1
	}else{
	  opserr << "WARNING invalid temptarget" << endln;
	  ifNoError = false;
	}

      	i += 1;


      } else if (strcmp(argv[i],"-physicalproperty")==0 && (i+8)<=(argc-1)) { //<-physicalproperty capap? capaf? capac? caparl? condp? condf? condc? condrl?> の読み込み
	
	if (Tcl_GetDouble(interp, argv[i+1], &capap) != TCL_OK || capap <= 0)  {
	  opserr << "WARNING invalid capap\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+2], &capaf) != TCL_OK || capaf <= 0)  {
	  opserr << "WARNING invalid capaf\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+3], &capac) != TCL_OK || capac <= 0)  {
	  opserr << "WARNING invalid capac\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+4], &caparl) != TCL_OK )  {
	  opserr << "WARNING invalid caparl\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+5], &condp) != TCL_OK || condp <= 0)  {
	  opserr << "WARNING invalid condp\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+6], &condf) != TCL_OK || condf <= 0)  {
	  opserr << "WARNING invalid condf\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+7], &condc) != TCL_OK || condc <= 0)  {
	  opserr << "WARNING invalid condc\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+8], &condrl) != TCL_OK )  {
	  opserr << "WARNING invalid condrl\n";
	  ifNoError = false;
	}

	i += 8;


      } else if (strcmp(argv[i],"-height")==0 && (i+6)<=(argc-1)) { // <-height hc? hf? hj? hrl? hss? hi?> の読み込み
	
	if (Tcl_GetDouble(interp, argv[i+1], &hc) != TCL_OK || hc <= 0)  {
	  opserr << "WARNING invalid hc\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+2], &hf) != TCL_OK || hf <= 0)  {
	  opserr << "WARNING invalid hf\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+3], &hj) != TCL_OK )  {
	  opserr << "WARNING invalid hj\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+4], &hrl) != TCL_OK )  {
	  opserr << "WARNING invalid hrl\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+5], &hss) != TCL_OK )  {
	  opserr << "WARNING invalid hss\n";
	  ifNoError = false;
	}

	if (Tcl_GetDouble(interp, argv[i+6], &hi) != TCL_OK )  {
	  opserr << "WARNING invalid hi\n";
	  ifNoError = false;
	}

	i += 6;



      } else if (strcmp(argv[i],"-loop")==0 && (i+1)<=(argc-1)) { // <-loop l?> の読み込み
	
	if (Tcl_GetInt(interp, argv[i+1], &loop) != TCL_OK || loop <= 0)  {
	  opserr << "WARNING invalid loop\n";
	  ifNoError = false;
	}

	i += 1;

      } else if (strcmp(argv[i],"-totalnumber")==0 && (i+1)<=(argc-1)) { // <-totalnumber n?> の読み込み
	
	if (Tcl_GetInt(interp, argv[i+1], &totalnumber) != TCL_OK || totalnumber <= 0)  {
	  opserr << "WARNING invalid totalnumber\n";
	  ifNoError = false;
	}

	i += 1;

      } else if (strcmp(argv[i],"-lim")==0 && (i+1)<=(argc-1)) { // <-lim limDisp?> の読み込み
	
	if (Tcl_GetDouble(interp, argv[i+1], &limDisp) != TCL_OK || limDisp < 0)  {
	  opserr << "WARNING invalid limDisp\n";
	  ifNoError = false;
	}

	i += 1;
	
      } else { //invalid option
	
	opserr << "WARNING invalid optional arguments \n";
	ifNoError = false;
	break;

      }
    }
  } //end input
  


  //confirm material
  if (recvMat != 1)  {
    opserr << "WARNING wrong number of -mat inputs\n";
    opserr << "got " << recvMat << " inputs, but want 1 input\n";
    ifNoError = false;
  }

  
  //if error detected
  if (!ifNoError) {
    //input:
    printCommand(argc, argv);
    //want:
    opserr << "Want: element MultipleShearSpring_forThermalConductivityAnalysis eleTag? iNode? jNode? nSpring? -mat matTag? <-lim dsp> <-dtstep dt?> <-mesh tate? yoko? hyoko? ctate? ftate? jtate? itate? norl? noss?> <-height hc? hf? hj? hrl? hss? hi?> <-temptarget targetorauto? targettate? targetyoko?> <-loop l?> <-totalnumber n?> <-physicalproperty capap? capaf? capac? caparl? condp? condf? condc? condrl?> <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>\n";
    return TCL_ERROR;
  }
  



  // now create the MultipleShearSpring_forThermalConductivityAnalysis
  theElement = new MultipleShearSpring_forThermalConductivityAnalysis(eleTag, iNode, jNode, nSpring, material, limDisp, dtstep, tate, yoko , hyoko, ctate, ftate, jtate, itate, norl, noss, hc, hf, hj, hrl, hss, hi, target_type, targettate, targetyoko, loop, totalnumber, capap, capaf, capac, caparl, condp, condf, condc, condrl, oriYp, oriX, mass);

  if (theElement == 0)  {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  // then add the MultipleShearSpring_forThermalConductivityAnalysis to the domain
  if (theTclDomain->addElement(theElement) == false)  {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }
  
  // if get here we have successfully created the MultipleShearSpring_forThermalConductivityAnalysis and added it to the domain
  return TCL_OK;
}


// initialize the class wide variables
Matrix MultipleShearSpring_forThermalConductivityAnalysis::theMatrix(12,12);
Vector MultipleShearSpring_forThermalConductivityAnalysis::theVector(12);
Vector MultipleShearSpring_forThermalConductivityAnalysis::theLoad(12);
Vector MultipleShearSpring_forThermalConductivityAnalysis::theHyst(7);




MultipleShearSpring_forThermalConductivityAnalysis::MultipleShearSpring_forThermalConductivityAnalysis(int Tag, int Nd1, int Nd2,
					 int NSpring,
				      KikuchiAikenLRB_TCA *Material,
                                      double LimDisp, double Dtstep, 
                                      int Tate, int Yoko, int Hyoko, int Ctate, int Ftate, int Jtate,  int Itate,
				      int Norl, int Noss,
				      double Hc, double Hf,  double Hj, double Hrl, double Hss, double Hi,
				      int Target_type,
				      int Targettate , int Targetyoko, int Loop, int Totalnumber,
                                      double Capap, double Capaf, double Capac, double Caparl,
	                              double Condp, double Condf, double Condc, double Condrl,
                                      const Vector OriYp, const Vector OriX, double Mass)
  : Element(Tag, ELE_TAG_MultipleShearSpring_forThermalConductivityAnalysis),
    connectedExternalNodes(2),
    nSpring(NSpring), limDisp(LimDisp), dtstep(Dtstep),
    tate(Tate), yoko(Yoko), hyoko(Hyoko), ctate(Ctate), ftate(Ftate), jtate(Jtate), itate(Itate),
    norl(Norl), noss(Noss),
    hc(Hc), hf(Hf), hj(Hj), hrl(Hrl), hss(Hss), hi(Hi),
    target_type(Target_type),
    targettate(Targettate), targetyoko(Targetyoko), loop(Loop), totalnumber(Totalnumber),
    capap(Capap), capaf(Capaf), capac(Capac), caparl(Caparl),
    condp(Condp), condf(Condf), condc(Condc), condrl(Condrl),
    allTemperatures(Tate*Yoko),
    oriX(OriX), oriYp(OriYp), mass(Mass),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6),
    basicQ2(6)
{
  

  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis::setUp() - element: "
	   << this->getTag() << " failed to create an ID of size 2\n";
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;
  
  // check material input
  if (Material == 0)  {
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis::MultipleShearSpring_forThermalConductivityAnalysis() - "
	   << "null uniaxial material pointer passed.\n";
    exit(-1);
  }

  theMaterials = new KikuchiAikenLRB_TCA* [nSpring];

  // get copies of the uniaxial materials
  for (int i=0; i<nSpring; i++)  {
    theMaterials[i] = (KikuchiAikenLRB_TCA *)Material->getCopy();

    if (theMaterials[i] == 0) {
      opserr << "MultipleShearSpring_forThermalConductivityAnalysis::MultipleShearSpring_forThermalConductivityAnalysis() - "
 	     << "failed to copy uniaxial material.\n";
      exit(-1);
    }
  }

  
  //arrangement of each spring
  cosTht = new double [nSpring];
  sinTht = new double [nSpring];
  
  for (int i=0; i<nSpring; i++)  {
    cosTht[i] = cos(M_PI*i/nSpring);
    sinTht[i] = sin(M_PI*i/nSpring);
  }

  trialF  = new double [2];
  commitF = new double [2];
  dDeform = new double [2];
  trialD  = new double [2];
  commitD = new double [2];


  TCA_data = new Data*[tate + 2];
  for(int i=0; i<tate+2; i++){
  TCA_data[i]  = new Data[yoko + 2];
  }


  //calculate Feq and Seq

  //imaginary material to caluculate Feq and Seq
  dmyMssMaterial = (KikuchiAikenLRB_TCA *)Material->getCopy();
  if (dmyMssMaterial == 0) {
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis::MultipleShearSpring_forThermalConductivityAnalysis() - "
	   << "failed to copy uniaxial material.\n";
    exit(-1);
  }
  dmyMssMaterial->revertToStart();


  //initial Feq and Seq
  if (limDisp > 0) {
    double uRef, fRef, sRef;//u:deformation, f:force, s:stiffness
    double uCmp, fSum, sSum;

    //imaginary material (1-directional deformation)
    uRef = limDisp;
    dmyMssMaterial->setTrialStrain(uRef,0);
    fRef = dmyMssMaterial->getStress();
    sRef = dmyMssMaterial->getTangent();

    //MSS
    fSum = 0.0;
    sSum = 0.0;
    for (int i=0; i<nSpring; i++)  {
      uCmp = uRef * cosTht[i];
      dmyMssMaterial->setTrialStrain(uCmp,0);
      fSum += dmyMssMaterial->getStress()  * cosTht[i];
      sSum += dmyMssMaterial->getTangent() * cosTht[i] * cosTht[i];
    }
    mssFeq = fRef/fSum;
    mssSeq = sRef/sSum;

  } else {

    mssFeq = 1.0;
    mssSeq = 1.0;

  }


  // initial basic stiffness matrix
  basicStiffInit.Zero();
  for (int i=0; i<nSpring; i++)  {
    double tmpTangent = theMaterials[i]->getInitialTangent();

    basicStiffInit(1,1) += tmpTangent * cosTht[i] * cosTht[i];
    basicStiffInit(1,2) += tmpTangent * cosTht[i] * sinTht[i];
    basicStiffInit(2,1) += tmpTangent * sinTht[i] * cosTht[i];
    basicStiffInit(2,2) += tmpTangent * sinTht[i] * sinTht[i];
  }


  // Seq: equivalent coefficient for stiffness
  basicStiffInit *= mssSeq;

  basicStiffInit *= totalnumber;//解析台数分の初期剛性を与える

  // initialize variables
  this->revertToStart();

}


MultipleShearSpring_forThermalConductivityAnalysis::MultipleShearSpring_forThermalConductivityAnalysis()
  : Element(0, ELE_TAG_MultipleShearSpring_forThermalConductivityAnalysis),
    connectedExternalNodes(2),
    nSpring(0), limDisp(0), dtstep(0.0),
    tate(0), yoko(0), hyoko(0), ctate(0), ftate(0), jtate(0), itate(0),
    norl(0), noss(0),
    hc(0.0), hf(0.0), hj(0.0), hrl(0.0), hss(0.0), hi(0.0),
    target_type(0),
    targettate(0), targetyoko(0), loop(0), totalnumber(0),
    capap(0.0), capaf(0.0), capac(0.0), caparl(0.0),
    condp(0.0), condf(0.0), condc(0.0), condrl(0.0),
    allTemperatures(0),
    oriX(0), oriYp(0), mass(0.0),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6),
    basicQ2(6)
{


  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis::MultipleShearSpring_forThermalConductivityAnalysis() - "
	   <<  "failed to create an ID of size 2\n";
    exit(-1);
  }
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  // set material array pointers to NULL
  theMaterials = 0;
  dmyMssMaterial = 0;

}


MultipleShearSpring_forThermalConductivityAnalysis::~MultipleShearSpring_forThermalConductivityAnalysis()
{
  // invoke the destructor on any objects created by the object
  // that the object still holds a pointer to

  if (theMaterials != 0) {
    for (int i=0; i<nSpring; i++)
      if (theMaterials[i] != 0)
	delete theMaterials[i];
    delete [] theMaterials;
  }

  if (cosTht != 0)
    delete [] cosTht;

  if (sinTht != 0)
    delete [] sinTht;

  if (dmyMssMaterial != 0)
    delete dmyMssMaterial;


  if (trialF != 0)
    delete [] trialF;

  if (commitF != 0)
    delete [] commitF;

  if (dDeform != 0)
    delete [] dDeform;

  if (trialD != 0)
    delete [] trialD;

  if (commitD != 0)
    delete [] commitD;

  for(int i=0; i<tate+2; i++){
    delete [] TCA_data[i];
  }
  delete TCA_data;



}


int MultipleShearSpring_forThermalConductivityAnalysis::getNumExternalNodes() const
{
  return 2;
}


const ID& MultipleShearSpring_forThermalConductivityAnalysis::getExternalNodes() 
{
  return connectedExternalNodes;
}


Node** MultipleShearSpring_forThermalConductivityAnalysis::getNodePtrs() 
{
  return theNodes;
}


int MultipleShearSpring_forThermalConductivityAnalysis::getNumDOF() 
{
  return 12;
}


void MultipleShearSpring_forThermalConductivityAnalysis::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (!theDomain)  {
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    return;
  }
  
  // first set the node pointers
  theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
  theNodes[1] = theDomain->getNode(connectedExternalNodes(1));	
  
  // if can't find both - send a warning message
  if (!theNodes[0] || !theNodes[1])  {
    if (!theNodes[0])  {
      opserr << "WARNING MultipleShearSpring_forThermalConductivityAnalysis::setDomain() - Nd1: " 
	     << connectedExternalNodes(0) << " does not exist in the model for ";
    } else  {
      opserr << "WARNING MultipleShearSpring_forThermalConductivityAnalysis::setDomain() - Nd2: " 
	     << connectedExternalNodes(1) << " does not exist in the model for ";
    }
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis ele: " << this->getTag() << endln;
    
    return;
  }
  
  // now determine the number of dof and the dimension    
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  
  // if differing dof at the ends - print a warning message
  if (dofNd1 != 6)  {
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis::setDomain() - node 1: "
	   << connectedExternalNodes(0) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  if (dofNd2 != 6)  {
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis::setDomain() - node 2: "
	   << connectedExternalNodes(1) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  
  // call the base class method
  this->DomainComponent::setDomain(theDomain);
  
  // set up the transformation matrix for orientation
  this->setUp();
}   	 


int MultipleShearSpring_forThermalConductivityAnalysis::commitState()
{

  int errCode = 0;
  
  // commit material models
  for (int i=0; i<nSpring; i++)
    errCode += theMaterials[i]->commitState();

  this->SumEnergy();// 履歴吸収エネルギー

  commitEsum =  trialEsum;

  theHyst(0) =  commitEsum;//commitEsumをはき出す

  this->ThermalConductivityAnalysis();//熱伝導解析


  for (int i=0; i<2; i++){
  commitF[i] =  trialF[i];
  commitD[i] =  trialD[i];
  }

  for (int i=0; i<nSpring; i++){
    errCode += theMaterials[i]->setTemp(target_temp);
  }

  dmyMssMaterial->setTemp(target_temp);

  return errCode;
}


int 
MultipleShearSpring_forThermalConductivityAnalysis::SumEnergy()
{
  delta_Esum = 0.0;

   for (int i=0; i<2; i++){
     trialD[i] = basicDisp(i+1);
     trialF[i] = basicQ2(i+1);
    dDeform[i] = trialD[i] - commitD[i];
    delta_Esum += ((trialF[i]+commitF[i])*dDeform[i] ) * 0.5;
  }

  trialEsum = commitEsum + delta_Esum;

  return 0;
}




int 
MultipleShearSpring_forThermalConductivityAnalysis::ThermalConductivityAnalysis()
{

 double dtstep = ops_Dt;
 dt = dtstep / loop;

    for(int n=1; n<=loop; n++){
      //計算部分
      for(int h=1; h<=tate; h++){
        for(int r=1; r<=yoko; r++){
          if(n == 1 && r<=hyoko && h>ctate+ftate+jtate && h<=ctate+ftate+jtate+rltate+sstate && ((h-(ctate+ftate+jtate)) % 2 == 1) ){
	    dfE = delta_Esum * 0.5 * TCA_data[h][r].volume / heated_part_volume;
          }else{
	    dfE = 0.0;
          }
          TCA_data[h][r].Trial_tempX = TCA_data[h][r].Commit_tempX
	    + dt/TCA_data[h][r].heat_capacity * (TCA_data[h-1][r].boundary_bottom * (TCA_data[h-1][r].Commit_tempX - TCA_data[h][r].Commit_tempX)
					- TCA_data[h][r].boundary_bottom  * (TCA_data[h][r].Commit_tempX   - TCA_data[h+1][r].Commit_tempX)
					+ TCA_data[h][r-1].boundary_right * (TCA_data[h][r-1].Commit_tempX - TCA_data[h][r].Commit_tempX)
					- TCA_data[h][r].boundary_right   * (TCA_data[h][r].Commit_tempX   - TCA_data[h][r+1].Commit_tempX))
	    + dfE / TCA_data[h][r].heat_capacity;


        //X値を代入する
	  TCA_data[h][r].Commit_tempX = TCA_data[h][r].Trial_tempX;
        }
      }
    }
    for(int h=1; h<=tate; h++){
      for(int r=1; r<=yoko; r++){
        TCA_data[h][r].Commit_temp = TCA_data[h][r].Commit_tempX;
        TCA_data[h][r].Trial_temp  = TCA_data[h][r].Trial_tempX;
        TCA_data[h][r].Commit_temp = TCA_data[h][r].Trial_temp;
      }
    }



    if(target_type == 0 ){//メッシュ固定

      target_temp = TCA_data[targettate][targetyoko].Trial_temp-273.0 ;

    }else if(target_type == 1){//最大値が採用温度

      double search_temp=0.0;


    // //(1)全鉛プラグから最大温度を探し出す場合
    // for(int h=ctate+ftate+1; h<=ctate+ftate+jtate+rltate+sstate+itate; h++){
    //   for(int r=1; r<=hyoko; r++){

    // 	search_temp = TCA_data[h][r].Trial_temp-273.0 ;

    // 	if(target_temp < search_temp )

    // 	  target_temp = search_temp;
    //   }
    // }

    //(2)鹿島の式を用いる場合
      target_temp = -273.0;
      for(int h=ctate+ftate+1; h<=ctate+ftate+jtate+rltate+sstate+itate; h++){

	search_temp = TCA_data[h][hyoko].Trial_temp-273.0 ;

	if(target_temp < search_temp )

	  target_temp = search_temp;

    }


    }else{

	opserr << "WARNING invalid target_type \n";
    }


    //------------温度をはき出す(仮に作成)-----------------------------
    theHyst(1) =  TCA_data[3][1].Trial_temp   -273.0 ;//鉛頂部温度
    theHyst(2) =  TCA_data[36][4].Trial_temp  -273.0 ;//鉛中心温度(詳細モデル用)
    theHyst(3) =  TCA_data[2][2].Trial_temp   -273.0 ;//フランジ中央温度
    theHyst(4) =  TCA_data[2][25].Trial_temp  -273.0 ;//フランジ端部温度
    theHyst(5) =  TCA_data[36][25].Trial_temp -273.0 ;//ゴム表面温度(詳細モデル用)
    theHyst(6) =  theMaterials[0]->getQd()  ;//降伏荷重
    //-----------------------------------------------------------------


    //------------全温度をはき出す-----------------------------
    for(int h=1; h<=tate; h++){
      for(int r=1; r<=yoko; r++){
        allTemperatures(yoko*(h-1)+(r-1)) = TCA_data[h][r].Trial_temp -273.0 ;//全温度をはき出す
      }
    }
    //------------全温度をはき出す-----------------------------


  return 0;
}


int MultipleShearSpring_forThermalConductivityAnalysis::revertToLastCommit()
{
  int errCode = 0;
  
  // revert material models
  for (int i=0; i<nSpring; i++)
    errCode += theMaterials[i]->revertToLastCommit();

  trialEsum =  commitEsum;

  for (int i=0; i<2; i++){
  trialF[i] =  commitF[i];
  trialD[i] =  commitD[i];
  }


  return errCode;
}


int MultipleShearSpring_forThermalConductivityAnalysis::revertToStart()
{
 //値の初期化(step0)

  int errCode=0;


  // trial variables
  basicDisp.Zero();
  basicForce.Zero();


  //---追加------------------------------------------------------------
  basicQ2.Zero();
  //-------------------------------------------------------------------


  // reset basic stiffness matrix
  basicStiff = basicStiffInit;
  
  // revert material models
  for (int i=0; i<nSpring; i++)
    errCode += theMaterials[i]->revertToStart();


  commitEsum = 0.0;
  trialEsum =  0.0;

  for (int i=0; i<2; i++){
  trialF[i]  = 0.0;
  commitF[i] = 0.0;
  trialD[i]  = 0.0;
  commitD[i] = 0.0;
  dDeform[i] = 0.0;
  }

  center_distance_right = 0.0;
  center_distance_bottom = 0.0;

  touch_area_right = 0.0;
  touch_area_bottom = 0.0;

  thermal_conductivity_right = 0.0;
  thermal_conductivity_bottom = 0.0;

  delta_Esum = 0.0;
  heated_part_ruiseki_r = 0.0;//鉛メッシュの半径情報の総和
  heated_part_ruiseki_h = 0.0;//鉛メッシュの高さ情報の総和
  heated_part_r = 0.0;//発熱部分（鉛）の半径
  heated_part_h = 0.0;//発熱部分（鉛）の総高さ
  heated_part_volume = 0.0;//発熱部分（鉛）の体積

  target_temp = 0.0;

  for(int h=0; h<=tate+1; h++){
    for(int r=0; r<=yoko+1; r++){
      TCA_data[h][r].Commit_temp  = theMaterials[0]->getTemp() + 273.0;  //初期温度はセ氏である
      TCA_data[h][r].Trial_temp   = theMaterials[0]->getTemp() + 273.0;    //値の初期化 
      TCA_data[h][r].Commit_tempX = theMaterials[0]->getTemp() + 273.0;  //初期温度はセ氏である
      TCA_data[h][r].Trial_tempX  = theMaterials[0]->getTemp() + 273.0;    //値の初期化 

      TCA_data[h][r].mesh_dr =1.0;
      TCA_data[h][r].mesh_dh =1.0;
      TCA_data[h][r].mesh_mat=0;
      TCA_data[h][r].heat_capacity=1.0;
      TCA_data[h][r].ruiseki_r=0.0;
      TCA_data[h][r].ruiseki_h=0.0;
      TCA_data[h][r].original_thermal_conductivity=0.0;
      TCA_data[h][r].volume=0.0;
      TCA_data[h][r].heat_capacity_by_volume=0.0;
      TCA_data[h][r].boundary_right=0.0;
      TCA_data[h][r].boundary_bottom=0.0;

    }
  }

  zentai_r = pow(theMaterials[0]->getAr() + theMaterials[0]->getAp() / M_PI , 0.5);//2016.03.10修正
  lead_r   = pow(theMaterials[0]->getAp() / M_PI , 0.5);


//天然ゴムと内部鋼板分離版

//諸元はここでファイルから読み込む(step1)
// [1][][] = meshのΔr[m], [2][][] = meshのΔh[m], [3][][] =物質識別番号(0= 鉛, 1= 天然ゴム, 2= フランジ or 挿入鋼板 or 連結鋼板 or 内部鋼板, 3=コンクリート)


//準備計算
    if(norl % 2 == 1){ //天然ゴム数が奇数の時
      rltate = norl / 2 + 1;
      sstate = noss / 2;
    }else{             //天然ゴム数が偶数の時
      rltate = norl / 2;
      sstate = noss / 2 + 1;
    }


  for(int h=1; h<=tate; h++){
    for(int r=1; r<=yoko; r++){


  //コンクリート部分
      if(h<=ctate){
        if(r<=hyoko){
	  TCA_data[h][r].mesh_dr=lead_r / hyoko;
        }else{
	  TCA_data[h][r].mesh_dr=(zentai_r - lead_r)/(yoko - hyoko);
        }
	TCA_data[h][r].mesh_dh=hc/ctate;
	TCA_data[h][r].mesh_mat=3;


 //フランジ部分
      }else if(h<=ctate+ftate){
        if(r<=hyoko){
	  TCA_data[h][r].mesh_dr=lead_r / hyoko;
        }else{
	  TCA_data[h][r].mesh_dr=(zentai_r - lead_r)/(yoko - hyoko);
        }
	TCA_data[h][r].mesh_dh=hf/ftate;
	TCA_data[h][r].mesh_mat=2;



 //連結鋼板部分
      }else if(h<=ctate+ftate+jtate && jtate>=1){
        if(r<=hyoko){
	  TCA_data[h][r].mesh_dr=lead_r / hyoko;
	  TCA_data[h][r].mesh_mat=0;
        }else{
	  TCA_data[h][r].mesh_dr=(zentai_r - lead_r)/(yoko - hyoko);
	  TCA_data[h][r].mesh_mat=2;
	}
	TCA_data[h][r].mesh_dh=hj/jtate;



 //積層ゴム＋隣接する鉛部分
      }else if(h<=ctate+ftate+jtate+rltate+sstate){

	//発熱鉛＋天然ゴム部分(積層ゴム奇数層目)
	if((h-(ctate+ftate+jtate)) % 2 == 1){

	  //発熱鉛部分
	  if(r<=hyoko){
	    TCA_data[h][r].mesh_dr=lead_r / hyoko;
	    TCA_data[h][r].mesh_mat=0;

	    if(h==ctate+ftate+jtate+rltate+sstate){//一番最後の要素は高さ半分
	      TCA_data[h][r].mesh_dh=hrl / 2;
	    }else{
	      TCA_data[h][r].mesh_dh=hrl;
	    }

	    //-------------------------------------------------------
	    heated_part_ruiseki_r =  heated_part_ruiseki_r + TCA_data[h][r].mesh_dr;
	    heated_part_ruiseki_h =  heated_part_ruiseki_h + TCA_data[h][r].mesh_dh;
	    //-------------------------------------------------------


	  }else{//天然ゴム部分
	    TCA_data[h][r].mesh_dr=(zentai_r - lead_r)/(yoko - hyoko);
	    TCA_data[h][r].mesh_mat=1;

	    if(h==ctate+ftate+jtate+rltate+sstate){//一番最後の要素は高さ半分
	      TCA_data[h][r].mesh_dh=hrl / 2;
	    }else{
	      TCA_data[h][r].mesh_dh=hrl;
	    }

	  }


	//非発熱鉛＋内部鋼板部分(積層ゴム偶数層目)
	}else if((h-(ctate+ftate+jtate)) % 2 == 0){

	  //非発熱鉛部分
	  if(r<=hyoko){
	    TCA_data[h][r].mesh_dr=lead_r / hyoko;
	    TCA_data[h][r].mesh_mat=0;

	    if(h==ctate+ftate+jtate+rltate+sstate){//一番最後の要素は高さ半分
	      TCA_data[h][r].mesh_dh=hss / 2;
	    }else{
	      TCA_data[h][r].mesh_dh=hss;
	    }


	  }else{//内部鋼板部分
	    TCA_data[h][r].mesh_dr=(zentai_r - lead_r)/(yoko - hyoko);
	    TCA_data[h][r].mesh_mat=2;

	    if(h==ctate+ftate+jtate+rltate+sstate){//一番最後の要素は高さ半分
	      TCA_data[h][r].mesh_dh=hss / 2;
	    }else{
	      TCA_data[h][r].mesh_dh=hss;
	    }

	  }
	}else{
	    opserr << "WARNING invalid step1 rl or ss" << endln;
	}



 //挿入鋼板部分
      }else if(h<=ctate+ftate+jtate+rltate+sstate+itate && itate>=1){
        if(r<=hyoko){
	  TCA_data[h][r].mesh_dr=lead_r / hyoko;
	  TCA_data[h][r].mesh_mat=0;

        }else{
	  TCA_data[h][r].mesh_dr=(zentai_r - lead_r)/(yoko - hyoko);
	  TCA_data[h][r].mesh_mat=2;
        }
	TCA_data[h][r].mesh_dh=hi*0.5/itate;


 //エラーの場合ここを通る
      }else{
	opserr << "WARNING invalid mesh step1" << endln;
      }



 //物性値をここで読み込む(step2)

    //鉛の物性値------------------
    if(TCA_data[h][r].mesh_mat == 0){
      TCA_data[h][r].heat_capacity_by_volume = capap;         //鉛の熱容量/体積[J/(m^3 * K)]
      TCA_data[h][r].original_thermal_conductivity = condp;   //鉛の熱伝導率[W/(m * K)]

    //天然ゴムの物性値-------
    }else if(TCA_data[h][r].mesh_mat == 1){
      TCA_data[h][r].heat_capacity_by_volume =  caparl;        //天然ゴムの熱容量/体積[J/(m^3 * K)]
      TCA_data[h][r].original_thermal_conductivity = condrl;   //天然ゴムの熱伝導率[W/(m * K)]

    //フランジ(SS400)の物性値-------
    }else if(TCA_data[h][r].mesh_mat == 2){
      TCA_data[h][r].heat_capacity_by_volume =  capaf;        //フランジの熱容量/体積[J/(m^3 * K)]
      TCA_data[h][r].original_thermal_conductivity = condf;   //フランジの熱伝導率[W/(m * K)]

    //コンクリートの物性値-------
    }else if(TCA_data[h][r].mesh_mat == 3){
      TCA_data[h][r].heat_capacity_by_volume =  capac;        //コンクリートの熱容量/体積[J/(m^3 * K)]
      TCA_data[h][r].original_thermal_conductivity = condc;   //コンクリートの熱伝導率[W/(m * K)]



    //エラーの場合ここを通る
    }else{
	opserr << "WARNING invalid mesh step2" << endln;
    }


 //heat_capacityの計算(step3)
    TCA_data[h][r].ruiseki_r = TCA_data[h][r-1].ruiseki_r+ TCA_data[h][r].mesh_dr;
    TCA_data[h][r].ruiseki_h = TCA_data[h-1][r].ruiseki_h+ TCA_data[h][r].mesh_dh;
    TCA_data[h][r].volume = (TCA_data[h][r].ruiseki_r * TCA_data[h][r].ruiseki_r - TCA_data[h][r-1].ruiseki_r * TCA_data[h][r-1].ruiseki_r)* M_PI * TCA_data[h][r].mesh_dh;
    TCA_data[h][r].heat_capacity = TCA_data[h][r].volume * TCA_data[h][r].heat_capacity_by_volume ;

    }
  }


 //発熱部分の体積の計算(step4)
  heated_part_r = heated_part_ruiseki_r / rltate;
  heated_part_h = heated_part_ruiseki_h / hyoko;
  heated_part_volume = heated_part_r *  heated_part_r * M_PI *  heated_part_h;


 //残りの諸計算(step5)
  for(int h=1; h<= tate; h++){
    for(int r=1; r<= yoko-1; r++){
      center_distance_right = (TCA_data[h][r].mesh_dr + TCA_data[h][r+1].mesh_dr) * 0.5;
      touch_area_right = 2.0 * M_PI * TCA_data[h][r].ruiseki_r * TCA_data[h][r].mesh_dh  ;
      thermal_conductivity_right = (TCA_data[h][r].original_thermal_conductivity + TCA_data[h][r+1].original_thermal_conductivity) * 0.5;
      TCA_data[h][r].boundary_right = touch_area_right * thermal_conductivity_right  / center_distance_right;
    }
  }
  for(int h=1; h<= tate-1; h++){
    for(int r=1; r<= yoko; r++){
      center_distance_bottom = (TCA_data[h][r].mesh_dh + TCA_data[h+1][r].mesh_dh) * 0.5;
      touch_area_bottom = (TCA_data[h][r].ruiseki_r * TCA_data[h][r].ruiseki_r - TCA_data[h][r-1].ruiseki_r * TCA_data[h][r-1].ruiseki_r)* M_PI ;
      thermal_conductivity_bottom = (TCA_data[h][r].original_thermal_conductivity + TCA_data[h+1][r].original_thermal_conductivity) * 0.5;
      TCA_data[h][r].boundary_bottom = touch_area_bottom * thermal_conductivity_bottom / center_distance_bottom;
    }
  }


  return errCode;
}


int MultipleShearSpring_forThermalConductivityAnalysis::update()
{


  // get global trial displacements and velocities
  const Vector &dsp1 = theNodes[0]->getTrialDisp();
  const Vector &dsp2 = theNodes[1]->getTrialDisp();
  const Vector &vel1 = theNodes[0]->getTrialVel();
  const Vector &vel2 = theNodes[1]->getTrialVel();
  
  static Vector globalDisp(12), globalDispDot(12);
  for (int i=0; i<6; i++)  {
    globalDisp(i)   = dsp1(i);  globalDispDot(i)   = vel1(i);
    globalDisp(i+6) = dsp2(i);  globalDispDot(i+6) = vel2(i);
  }

  static Vector localDispDot(12), basicDispDot(6);

  
  // transform response from the global to the local system
  localDisp    = Tgl*globalDisp;
  localDispDot = Tgl*globalDispDot;
  
  // transform response from the local to the basic system
  basicDisp    = Tlb*localDisp;
  basicDispDot = Tlb*localDispDot;


  // calculate shear forces and stiffnesses in basic y- and z-direction
  // get trial shear forces of hysteretic component
  basicForce.Zero();
  basicStiff.Zero();


  //---追加------------------------------------------------------------
  basicQ2.Zero();
  //-------------------------------------------------------------------


  for (int i=0; i<nSpring; i++)  {
    double tmpStrain     = basicDisp(1)*cosTht[i]    + basicDisp(2)*sinTht[i];
    double tmpStrainRate = basicDispDot(1)*cosTht[i] + basicDispDot(2)*sinTht[i];

    theMaterials[i]->setTrialStrain(tmpStrain,tmpStrainRate);
    double tmpStress  = theMaterials[i]->getStress();


    basicForce(1) += tmpStress * cosTht[i]; //basic-y
    basicForce(2) += tmpStress * sinTht[i]; //basic-z


  //---追加------------------------------------------------------------
    double tmpQ2  = theMaterials[i]->getQ2();

    basicQ2(1) += tmpQ2 * cosTht[i]; //basic-y
    basicQ2(2) += tmpQ2 * sinTht[i]; //basic-z
  //-------------------------------------------------------------------

    double tmpTangent = theMaterials[i]->getTangent();

    basicStiff(1,1) += tmpTangent * cosTht[i] * cosTht[i];
    basicStiff(1,2) += tmpTangent * cosTht[i] * sinTht[i];
    basicStiff(2,1) += tmpTangent * sinTht[i] * cosTht[i];
    basicStiff(2,2) += tmpTangent * sinTht[i] * sinTht[i];
  }


  // calculate Feq and Seq
  if (limDisp > 0) {
    double uRef, fRef, sRef;//u:deformation, f:force, s:stiffness
    double uCmp, fSum, sSum;
    double refDisp;

    //imaginary material (1-directional deformation)
    refDisp = sqrt(basicDisp(1)*basicDisp(1)+basicDisp(2)*basicDisp(2));
    uRef = (refDisp>limDisp) ? refDisp : limDisp;

    dmyMssMaterial->setTrialStrain(uRef,0);
    fRef = dmyMssMaterial->getStress();
    sRef = dmyMssMaterial->getTangent();

    //MSS
    fSum = 0.0;
    sSum = 0.0;
    for (int i=0; i<nSpring; i++)  {
      uCmp = uRef * cosTht[i];
      dmyMssMaterial->setTrialStrain(uCmp,0);
      fSum += dmyMssMaterial->getStress() * cosTht[i];
      sSum += dmyMssMaterial->getTangent() * cosTht[i] * cosTht[i];
    }
    mssFeq = fRef/fSum;
    mssSeq = sRef/sSum;

  }

  basicForce *= mssFeq;
  basicStiff *= mssSeq;


  basicQ2 *= mssFeq;//追加

  basicForce *= totalnumber;//解析台数分の荷重と剛性を返す
  basicStiff *= totalnumber;//解析台数分の荷重と剛性を返す

  basicQ2 *= totalnumber;//LRBの員数オプション用(解析の単純化)


  return 0;
}


double MultipleShearSpring_forThermalConductivityAnalysis::temp(void)
{
  return target_temp;
}


const Matrix& MultipleShearSpring_forThermalConductivityAnalysis::getTangentStiff()
{
  // zero the matrix
  theMatrix.Zero();
  
  // transform from basic to local system
  static Matrix localStiff(12,12);
  localStiff.addMatrixTripleProduct(0.0, Tlb, basicStiff, 1.0);
  
  // transform from local to global system
  theMatrix.addMatrixTripleProduct(0.0, Tgl, localStiff, 1.0);
  
  return theMatrix;
}


const Matrix& MultipleShearSpring_forThermalConductivityAnalysis::getInitialStiff()
{
  // zero the matrix
  theMatrix.Zero();
  
  // transform from basic to local system
  static Matrix localStiff(12,12);
  localStiff.addMatrixTripleProduct(0.0, Tlb, basicStiffInit, 1.0);
  
  // transform from local to global system
  theMatrix.addMatrixTripleProduct(0.0, Tgl, localStiff, 1.0);
  
  return theMatrix;
}


const Matrix& MultipleShearSpring_forThermalConductivityAnalysis::getMass()
{
  // zero the matrix
  theMatrix.Zero();
  
  // check for quick return
  if (mass == 0.0)  {
    return theMatrix;
  }    
  
  double m = 0.5*mass;
  for (int i = 0; i < 3; i++)  {
    theMatrix(i,i)     = m;
    theMatrix(i+6,i+6) = m;
  }
  
  return theMatrix; 
}


void MultipleShearSpring_forThermalConductivityAnalysis::zeroLoad()
{
  theLoad.Zero();
}


int MultipleShearSpring_forThermalConductivityAnalysis::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
  opserr <<"MultipleShearSpring_forThermalConductivityAnalysis::addLoad() - "
	 << "load type unknown for element: "
	 << this->getTag() << endln;
  
  return -1;
}


int MultipleShearSpring_forThermalConductivityAnalysis::addInertiaLoadToUnbalance(const Vector &accel)
{
  // check for quick return
  if (mass == 0.0)  {
    return 0;
  }
  
  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis::addInertiaLoadToUnbalance() - "
	   << "matrix and vector sizes are incompatible\n";
    return -1;
  }
  
  // want to add ( - fact * M R * accel ) to unbalance
  // take advantage of lumped mass matrix
  double m = 0.5*mass;
  for (int i = 0; i < 3; i++)  {
    theLoad(i)   -= m * Raccel1(i);
    theLoad(i+6) -= m * Raccel2(i);
  }
  
  return 0;
}


const Vector& MultipleShearSpring_forThermalConductivityAnalysis::getResistingForce()
{
  // zero the residual
  theVector.Zero();
  
  // determine resisting forces in local system
  static Vector localForce(12);
  localForce = Tlb^basicForce;
  
  // determine resisting forces in global system
  theVector = Tgl^localForce;
  
  // subtract external load
  theVector.addVector(1.0, theLoad, -1.0);
  
  return theVector;
}


const Vector& MultipleShearSpring_forThermalConductivityAnalysis::getResistingForceIncInertia()
{	
  theVector = this->getResistingForce();
  
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    theVector += this->getRayleighDampingForces();
  
  // now include the mass portion
  if (mass != 0.0)  {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();    
    
    double m = 0.5*mass;
    for (int i = 0; i < 3; i++)  {
      theVector(i)   += m * accel1(i);
      theVector(i+6) += m * accel2(i);
    }
  }
  
  return theVector;
}


int MultipleShearSpring_forThermalConductivityAnalysis::sendSelf(int commitTag, Channel &sChannel)
{
  return -1;
}


int MultipleShearSpring_forThermalConductivityAnalysis::recvSelf(int commitTag, Channel &rChannel,
				  FEM_ObjectBroker &theBroker)
{
  return -1;
}


int MultipleShearSpring_forThermalConductivityAnalysis::displaySelf(Renderer &theViewer,
				     int displayMode, float fact)
{
  // first determine the end points of the element based on
  // the display factor (a measure of the distorted image)
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	
  
  const Vector &end1Disp = theNodes[0]->getDisp();
  const Vector &end2Disp = theNodes[1]->getDisp();
  
  static Vector v1(3);
  static Vector v2(3);
  
  for (int i = 0; i < 3; i++)  {
    v1(i) = end1Crd(i) + end1Disp(i)*fact;
    v2(i) = end2Crd(i) + end2Disp(i)*fact;    
  }
  
  return theViewer.drawLine (v1, v2, 1.0, 1.0);
}


void MultipleShearSpring_forThermalConductivityAnalysis::Print(OPS_Stream &s, int flag)
{
  
  if (flag == 0)  {
    // print everything
    s << "Element: " << this->getTag(); 
    s << "  type: MultipleShearSpring_forThermalConductivityAnalysis  iNode: " << connectedExternalNodes(0);
    s << "  jNode: " << connectedExternalNodes(1) << endln;
    s << "  Material : " << theMaterials[0]->getTag() << endln;
    s << "  mass: " << mass << endln;
    // determine resisting forces in global system
    s << "  resisting force: " << this->getResistingForce() << endln;
  } else if (flag == 1)  {
    // does nothing
  }
}


Response* MultipleShearSpring_forThermalConductivityAnalysis::setResponse(const char **argv, int argc,
					   OPS_Stream &output)
{
  Response *theResponse = 0;
  
  output.tag("ElementOutput");
  output.attr("eleType","MultipleShearSpring_forThermalConductivityAnalysis");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  
  // global forces
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
    {
      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");
      
      theResponse = new ElementResponse(this, 1, theVector);
    }
  // local forces
  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
    {
      output.tag("ResponseType","N_ 1");
      output.tag("ResponseType","Vy_1");
      output.tag("ResponseType","Vz_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Tz_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","T_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");
      
      theResponse = new ElementResponse(this, 2, theVector);
    }
  // basic forces
  else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0)
    {
      output.tag("ResponseType","qb1");
      output.tag("ResponseType","qb2");
      output.tag("ResponseType","qb3");
      output.tag("ResponseType","qb4");
      output.tag("ResponseType","qb5");
      output.tag("ResponseType","qb6");
      
      theResponse = new ElementResponse(this, 3, Vector(6));
    }
  // local displacements
  else if (strcmp(argv[0],"localDisplacement") == 0 ||
	   strcmp(argv[0],"localDisplacements") == 0)
    {
      output.tag("ResponseType","ux_1");
      output.tag("ResponseType","uy_1");
      output.tag("ResponseType","uz_1");
      output.tag("ResponseType","rx_1");
      output.tag("ResponseType","ry_1");
      output.tag("ResponseType","rz_1");
      output.tag("ResponseType","ux_2");
      output.tag("ResponseType","uy_2");
      output.tag("ResponseType","uz_2");
      output.tag("ResponseType","rx_2");
      output.tag("ResponseType","ry_2");
      output.tag("ResponseType","rz_2");
      
      theResponse = new ElementResponse(this, 4, theVector);
    }
  // basic displacements
  else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"deformations") == 0 || 
	   strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0 ||
	   strcmp(argv[0],"basicDisplacement") == 0 || strcmp(argv[0],"basicDisplacements") == 0)
    {
      output.tag("ResponseType","ub1");
      output.tag("ResponseType","ub2");
      output.tag("ResponseType","ub3");
      output.tag("ResponseType","ub4");
      output.tag("ResponseType","ub5");
      output.tag("ResponseType","ub6");
      
      theResponse = new ElementResponse(this, 5, Vector(6));
    }
  

  // hysteresis info
  else if (strcmp(argv[0],"hyst") == 0 ||
           strcmp(argv[0],"hysteresis") == 0)
    {

      output.tag("ResponseType","commitEsum");
      output.tag("ResponseType","temp_lead1");
      output.tag("ResponseType","temp_lead2");
      output.tag("ResponseType","temp_steel1");
      output.tag("ResponseType","temp_steel2");
      output.tag("ResponseType","temp_rubber");
      output.tag("ResponseType","Qd");

      theResponse = new ElementResponse(this, 6, theHyst);
    }

  // allTemperatures info
  else if (strcmp(argv[0],"allTemperatures") == 0 ||
           strcmp(argv[0],"allTemp") == 0)
    {



    for(int h=1; h<=tate; h++){
      for(int r=1; r<=yoko; r++){
        char filename[1000];//仮の値
        sprintf(filename,"Temp_%d_%d",h,r);
        output.tag("ResponseType", filename );
      }
    }



      theResponse = new ElementResponse(this, 7, allTemperatures);
    }


  output.endTag(); // ElementOutput
  
  return theResponse;
}


int MultipleShearSpring_forThermalConductivityAnalysis::getResponse(int responseID, Information &eleInfo)
{
  switch (responseID)  {
  case 1:  // global forces
    return eleInfo.setVector(this->getResistingForce());
    
  case 2:  // local forces
    theVector.Zero();
    // determine resisting forces in local system
    theVector = Tlb^basicForce;
    
    return eleInfo.setVector(theVector);
    
  case 3:  // basic forces
    return eleInfo.setVector(basicForce);
    
  case 4:  // local displacements
    return eleInfo.setVector(localDisp);
    
  case 5:  // basic displacements
    return eleInfo.setVector(basicDisp);

  case 6:  // hysteresis information
    return eleInfo.setVector(theHyst);

  case 7:  // allTemperatures information
    return eleInfo.setVector(allTemperatures);

    
  default:
    return -1;
  }
}


// set up the transformation matrix for orientation
void MultipleShearSpring_forThermalConductivityAnalysis::setUp()
{ 
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	
  Vector oriXp = end2Crd - end1Crd;
  double elmLen = oriXp.Norm();
  
  if (elmLen > DBL_EPSILON)  {
    if (oriX.Size() == 0)  {
      oriX.resize(3);
      oriX = oriXp;
    } else  {
      opserr << "WARNING MultipleShearSpring_forThermalConductivityAnalysis::setUp() - " 
	     << "element: " << this->getTag() << endln
	     << "ignoring nodes and using specified "
	     << "local x vector to determine orientation\n";
    }
  }
  // check that vectors for orientation are of correct size
  if (oriX.Size() != 3 || oriYp.Size() != 3)  {
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis::setUp() - "
	   << "element: " << this->getTag() << endln
	   << "incorrect dimension of orientation vectors\n";
    exit(-1);
  }
    
  // establish orientation of element for the tranformation matrix
  // z = x cross yp
  Vector oriZ(3);
  oriZ(0) = oriX(1)*oriYp(2) - oriX(2)*oriYp(1);
  oriZ(1) = oriX(2)*oriYp(0) - oriX(0)*oriYp(2);
  oriZ(2) = oriX(0)*oriYp(1) - oriX(1)*oriYp(0);
  
  // y = z cross x
  Vector oriY(3);
  oriY(0) = oriZ(1)*oriX(2) - oriZ(2)*oriX(1);
  oriY(1) = oriZ(2)*oriX(0) - oriZ(0)*oriX(2);
  oriY(2) = oriZ(0)*oriX(1) - oriZ(1)*oriX(0);
    
  // compute length(norm) of vectors
  double xn = oriX.Norm();
  double yn = oriY.Norm();
  double zn = oriZ.Norm();
  
  // check valid x and y vectors, i.e. not parallel and of zero length
  if (xn == 0 || yn == 0 || zn == 0)  {
    opserr << "MultipleShearSpring_forThermalConductivityAnalysis::setUp() - "
	   << "element: " << this->getTag() << endln
	   << "invalid orientation vectors\n";
    exit(-1);
  }
  
  // create transformation matrix from global to local system
  Tgl.Zero();
  Tgl(0,0) = Tgl(3,3) = Tgl(6,6) = Tgl(9,9)   = oriX(0)/xn;
  Tgl(0,1) = Tgl(3,4) = Tgl(6,7) = Tgl(9,10)  = oriX(1)/xn;
  Tgl(0,2) = Tgl(3,5) = Tgl(6,8) = Tgl(9,11)  = oriX(2)/xn;
  Tgl(1,0) = Tgl(4,3) = Tgl(7,6) = Tgl(10,9)  = oriY(0)/yn;
  Tgl(1,1) = Tgl(4,4) = Tgl(7,7) = Tgl(10,10) = oriY(1)/yn;
  Tgl(1,2) = Tgl(4,5) = Tgl(7,8) = Tgl(10,11) = oriY(2)/yn;
  Tgl(2,0) = Tgl(5,3) = Tgl(8,6) = Tgl(11,9)  = oriZ(0)/zn;
  Tgl(2,1) = Tgl(5,4) = Tgl(8,7) = Tgl(11,10) = oriZ(1)/zn;
  Tgl(2,2) = Tgl(5,5) = Tgl(8,8) = Tgl(11,11) = oriZ(2)/zn;
  
  // create transformation matrix from local to basic system (linear)
  Tlb.Zero();
  Tlb(0,0) = Tlb(1,1) = Tlb(2,2) = Tlb(3,3) = Tlb(4,4) = Tlb(5,5) = -1.0;
  Tlb(0,6) = Tlb(1,7) = Tlb(2,8) = Tlb(3,9) = Tlb(4,10) = Tlb(5,11) = 1.0;
  Tlb(1,5) = Tlb(1,11) = -0.5*elmLen;
  Tlb(2,4) = Tlb(2,10) = 0.5*elmLen;
}

