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

// $Revision: 1.0 $
// $Date: 2014-02-17 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/EBSBearing/EBSBearing.h,v $
// Written: Ippei Tsunezawa
// Created: February 2014
// Modified: Feb 09, 2015
//
// Description: This file contains the implementation of the EBSBearing class.
//
// YamamotoBiaxialHDR.cppをベースにして作成
//
// local-x,local-y,local-z方向にばねを設定する
// その他の方向は剛性を0とする
// TclModelBuilderをコンストラクタの前に入れる。
//

#include <EBSBearing.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>

#include <ID.h>
#include <Vector.h>
#include <TclModelBuilder.h>

#include <float.h>
#define _USE_MATH_DEFINES  // M_PIをPIとして使えるようにする。この位置に置くこと。
#include <math.h>
#include <stdlib.h>
#include <string.h>


// --- TclModelBuilder_addEBSBearing の挿入 初め -------------------------------------------------

extern void printCommand(int argc, TCL_Char **argv);


int TclModelBuilder_addEBSBearing(ClientData clientData,
				      Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
				      TclModelBuilder *theTclBuilder)
{

  // ensure the destructor has not been called
  if (theTclBuilder == 0)  {
    opserr << "WARNING builder has been destroyed - EBSBearing\n";    
    return TCL_ERROR;
  }

  //自由度のチェック
  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING EBSBearing command only works when ndm is 3 and ndf is 6" << endln;
    return TCL_ERROR;
  }

  //最低限必要な引数
  int eleTag;
  int iNode;
  int jNode;

  int Tp;
  double DDo;
  double DDi;
  double Hr;
  double Mu;
  double Ec;
  double Gr;
  double Fz0;
  double As;

  //オプション引数のデフォルト値
  Vector oriX(0);
  Vector oriYp(3); oriYp(0) = 0.0; oriYp(1) = 1.0; oriYp(2) = 0.0;
  double mass = 0.0;


  //作成するオブジェクト
  Element *theElement = 0;

  //入力引数のエラーチェック
  bool ifNoError = true;


  if (argc < 14)  { //element EBSBearing eleTag? iNode? jNode? Tp? DDo? DDi? Hr? Mu? Ec? Gr? Fz0? As?
    // argc =            1           2             3      4      5     6   7    8    9    10     11       12       13        14
    // argv =       argv[0]      argv[1]      argv[2]  argv[3] ................. argv[8] argv[9] argv[10] argv[11] argv[12] argv[13]
    opserr << "WARNING insufficient arguments\n";
    ifNoError = false;

  } else {


    //argv[2~8]のチェック
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK)  {
      opserr << "WARNING invalid EBSBearing eleTag\n";
      ifNoError = false;
    }

    // iNode（i節点番号）の読み込み
    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK)  {
      opserr << "WARNING invalid iNode\n";
      ifNoError = false;
    }

    // jNode（j節点番号）の読み込み
    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK)  {
      opserr << "WARNING invalid jNode\n";
      ifNoError = false;
    }

    // Tp（積層ゴム復元力モデルのタイプ）の読み込み
    if (strcmp(argv[5],"1") == 0) {
      Tp = 1; //弾性すべり支承
    } else if (strcmp(argv[5],"2") == 0) {
      Tp = 2; // other
    } else {
      opserr << "WARNING invalid EBSBearing Tp" << endln;
      ifNoError = false;
    }

    // DDo（積層ゴム外径）の読み込み
    if (Tcl_GetDouble(interp, argv[6], &DDo) != TCL_OK || DDo <= 0.0) {
      opserr << "WARNING invalid EBSBearing DDo" << endln;
      ifNoError = false;
    }

    // DDi（積層ゴム内径）の読み込み
    if (Tcl_GetDouble(interp, argv[7], &DDi) != TCL_OK || DDi < 0.0) {
      opserr << "WARNING invalid EBSBearing DDi" << endln;
      ifNoError = false;
    }

    // Hr（ゴム総厚）の読み込み
    if (Tcl_GetDouble(interp, argv[8], &Hr) != TCL_OK || Hr <= 0.0) {
      opserr << "WARNING invalid EBSBearing Hr" << endln;
      ifNoError = false;
    }

    // Mu（摩擦係数）の読み込み
    if (Tcl_GetDouble(interp, argv[9], &Mu) != TCL_OK || Mu < 0.0) {
      opserr << "WARNING invalid EBSBearing Mu" << endln;
      ifNoError = false;
    }

    // Ec（圧縮弾性係数）の読み込み
    if (Tcl_GetDouble(interp, argv[10], &Ec) != TCL_OK || Ec < 0.0) {
      opserr << "WARNING invalid EBSBearing Ec" << endln;
      ifNoError = false;
    }

    // Gr（せん断弾性係数）の読み込み
    if (Tcl_GetDouble(interp, argv[11], &Gr) != TCL_OK || Gr < 0.0) {
      opserr << "WARNING invalid EBSBearing Gr" << endln;
      ifNoError = false;
    }

    // Fz0（初期軸力）の読み込み
    if (Tcl_GetDouble(interp, argv[12], &Fz0) != TCL_OK ) {
      opserr << "WARNING invalid EBSBearing Fz0" << endln;
      ifNoError = false;
    }

    // As（すべり材の断面積）の読み込み
    if (Tcl_GetDouble(interp, argv[13], &As) != TCL_OK ) {
      opserr << "WARNING invalid EBSBearing As" << endln;
      ifNoError = false;
    }


    // check print--------------------------------------------/
    //  opserr << "   \n";
    //  opserr << "TclModelBuilder_addEBSBearing()\n";
    //  opserr << "  tp  = " << Tp << endln;
    //  opserr << "  ddo = " << DDo << endln;
    //  opserr << "  ddi = " << DDi << endln;
    //  opserr << "  hr  = " << Hr << endln;
    //  opserr << "  mu  = " << Mu << endln;
    //  opserr << "  ec  = " << Ec << endln;
    //  opserr << "  gr  = " << Gr << endln;
    //  opserr << "  Fz0 = " << Fz0 << endln;
    //  opserr << "  As  = " << As  << endln;
    //------------------------------------------------------

    // argv[10~]のチェック
    for (int i=14; i<=(argc-1); i++) {
      double value;

      if (strcmp(argv[i],"-orient")==0 && (i+6)<=(argc-1) && Tcl_GetDouble(interp,argv[i+4], &value) == TCL_OK) { // <-orient x1? x2? x3? yp1? yp2? yp3?> の読み込み
	//                                                                     argv[i+4]が実数であることを判断                            yp1=argv[i+4]
	oriX.resize(3);  //  入力データが x1? x2? x3? yp1? yp2? yp3?のパターンなので、ori.X(3)の領域を確保

	// x1, x2, x3を読み込む
	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    opserr << "WARNING invalid -orient value\n";
	    ifNoError = false;
	  } else {
	    oriX(j-1) = value;
	  }
	}
	
	i += 3;
	
	// yp1, yp2, yp3を読み込む
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
	
	// yp1, yp2, yp3を読み込む
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
	
	// massを読み込む
	if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK || mass <= 0)  {
	  opserr << "WARNING invalid mass\n";
	  ifNoError = false;
	}
	
      } else { //無効なオプション
	
	opserr << "WARNING invalid optional arguments \n";
	ifNoError = false;
	break;
	
      }
    }

  } //引数の読み込み終了

  
  //入力にエラーがあった場合の処理
  if (!ifNoError) {
    //入力データ
    printCommand(argc, argv);
    //必要データ
    opserr << "Want: element EBSBearing eleTag? iNode? jNode? Tp? DDo? DDi? Hr? Mu? Ec? Gr? Fz0? As?<-orient <x1? x2? x3?> y1? y2? y3?> <-mass m?>\n";
    return TCL_ERROR;
  }
  

  // now create the EBSBearing
  theElement = new EBSBearing(eleTag, iNode, jNode, Tp, DDo, DDi, Hr, Mu, Ec, Gr, Fz0, As, oriYp, oriX, mass);

  if (theElement == 0)  {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "EBSBearing element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  // then add the EBSBearing to the domain
  if (theTclDomain->addElement(theElement) == false)  {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "EBSBearing element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }       
  
  // if get here we have successfully created the EBSBearing and added it to the domain
  return TCL_OK;
  }


// --- TclModelBuilder_addEBSBearing の挿入 終わり -----------------------------------------------


// initialize the class wide variables
Matrix EBSBearing::theMatrix(12,12);
Vector EBSBearing::theVector(12);
Vector EBSBearing::theLoad(12);


EBSBearing::EBSBearing(int Tag, int Nd1, int Nd2, int Tp, double DDo, double DDi, double Hr, double Mu, double Ec, double Gr, double Fz0, double As,
	   const Vector OriYp, const Vector OriX, double Mass)
  : Element(Tag, ELE_TAG_EBSBearing),
    tp(Tp),ddo(DDo),ddi(DDi),hr(Hr),mu(Mu),ec(Ec),gr(Gr),fz0(Fz0), as(As),
    connectedExternalNodes(2),
    oriX(OriX), oriYp(OriYp), mass(Mass),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6)
{
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "EBSBearing::setUp() - element: "
	   << this->getTag() << " failed to create an ID of size 2\n";
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;



  // 積層ゴムの断面積
  ar = M_PI * ( ddo*ddo - ddi*ddi ) / 4;



  // ゴム種ごとのパラメータ設定
  //弾性すべり支承
  if ( tp == 1 ) {

      initialStiff[0] = gr * ar / hr;
      initialStiff[1] = gr * ar / hr;
      initialStiff[2] = ec * ar / hr;
  }
  // // 弾性すべり支承(摩擦係数の面圧依存を考慮)
  if ( tp == 2 ) {

      initialStiff[0] = gr * ar / hr;
      initialStiff[1] = gr * ar / hr;
      initialStiff[2] = ec * ar / hr;
  }


  // initial basic stiffness matrix
  basicStiffInit.Zero();

  basicStiffInit(0,0) = this->getInitialTangent(2);
  basicStiffInit(1,1) = this->getInitialTangent(0);
  basicStiffInit(2,2) = this->getInitialTangent(1);



  // initialize variables
  this->revertToStart();

  //check print
  opserr << "basicStiffInit:  " << basicStiff << endln;


  // check print--------------------------------------------
  //opserr << "   \n";
  //opserr << "EBSBearing::EBSBearing()\n";
  //opserr << "  tp = " << tp << endln;
  //opserr << "  dr = " << dr << endln;
  //opserr << "  hr = " << hr << endln;
  // opserr << "  ar = " << ar << endln;
  // opserr << "  M_PI = " << M_PI << endln;
  // opserr << "  ddo = " << ddo << endln;
  // opserr << "  ddi = " << ddi << endln;
  //opserr << "  nn = " << nn << endln;
  //------------------------------------------------------

}

EBSBearing::EBSBearing()
  : Element(0, ELE_TAG_EBSBearing),
    connectedExternalNodes(2),
    oriX(0), oriYp(0), mass(0.0),
    Tgl(12,12), Tlb(6,12),
    basicDisp(6), localDisp(12), basicForce(6), basicStiff(6,6), basicStiffInit(6,6)
{	

  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "EBSBearing::EBSBearing() - "
	   <<  "failed to create an ID of size 2\n";
    exit(-1);
  }
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  //
  for (int i=0; i<3; i++) {

    trialDeform[i]  = 0.0;
    trialForce[i]   = 0.0;
    trialStiff[i]   = 0.0;
    trialQ[i]   = 0.0;
//    trialP[i]   = 0.0;

    commitDeform[i]  = 0.0;
    commitForce[i]   = 0.0;
    commitStiff[i]   = 0.0;
    commitQ[i]   = 0.0;

  }
}



EBSBearing::~EBSBearing()
{
  // does nothing
}


int EBSBearing::getNumExternalNodes() const
{
  return 2;
}


const ID& EBSBearing::getExternalNodes() 
{
  return connectedExternalNodes;
}


Node** EBSBearing::getNodePtrs() 
{
  return theNodes;
}


int EBSBearing::getNumDOF() 
{
  return 12;
}


void EBSBearing::setDomain(Domain *theDomain)
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
      opserr << "WARNING EBSBearing::setDomain() - Nd1: " 
	     << connectedExternalNodes(0) << " does not exist in the model for ";
    } else  {
      opserr << "WARNING EBSBearing::setDomain() - Nd2: " 
	     << connectedExternalNodes(1) << " does not exist in the model for ";
    }
    opserr << "EBSBearing ele: " << this->getTag() << endln;
    
    return;
  }
  
  // now determine the number of dof and the dimension    
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  
  // if differing dof at the ends - print a warning message
  if (dofNd1 != 6)  {
    opserr << "EBSBearing::setDomain() - node 1: "
	   << connectedExternalNodes(0) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  if (dofNd2 != 6)  {
    opserr << "EBSBearing::setDomain() - node 2: "
	   << connectedExternalNodes(1) << " has incorrect number of DOF (not 6)\n";
    return;
  }
  
  // call the base class method
  this->DomainComponent::setDomain(theDomain);
  
  // set up the transformation matrix for orientation
  this->setUp();
}   	 


int EBSBearing::setTrialStrain(const Vector &strain)
{

  trialDeform[0] = strain(1); // local-y
  trialDeform[1] = strain(2); // local-z
  trialDeform[2] = strain(0); // local-x,2014.02.24 EBSBearingの挿入
      
  Kz = ec * ar / hr;//Kzの算出
  Kxy = gr * ar / hr;//Kxyの算出

  // check print--------------------------------------------
  //opserr << "   \n";
  //opserr << "EBSBearing::setTrialStrain()\n";
  //opserr << "  trialDeform[0] = " << trialDeform[0] << endln;
  //opserr << "  trialDeform[1] = " << trialDeform[1] << endln;
  //opserr << "  trialP[0] = " << trialP[0] << endln;
  //opserr << "  trialP[1] = " << trialP[1] << endln;
  //------------------------------------------------------


  if ( tp == 1 ) {
  //-----------------------------------------------------------------------------
  //弾性すべり支承
  //-----------------------------------------------------------------------------

      if(trialDeform[2]  >= -fz0 / Kz){
      trialForce[2] = -fz0;
      Alpha = 0.0;//Alphaの算出
      }else{
      trialForce[2] = Kz * trialDeform[2];
      Alpha = mu * (-trialForce[2] -fz0) / Kxy;//Alphaの算出

      }
  }  //   if ( tp == 1 ) ... end-----------------------------------------------


  else if ( tp == 2 ) {
  //-----------------------------------------------------------------------------
  //弾性すべり支承(摩擦係数の面圧依存を考慮)
  //-----------------------------------------------------------------------------

      if(trialDeform[2]  >= -fz0 / Kz){      
      trialForce[2] = -fz0;
      mui = 0.0;//muiの算出
      Alpha = 0.0;//Alphaの算出
      }else{
      trialForce[2] = Kz * trialDeform[2];
      mui = 0.305 * pow((-trialForce[2] -fz0) / as * 1.0e-6 ,-0.29); //muiの算出
      Alpha = mui * (-trialForce[2] -fz0) / Kxy;//Alphaの算出
      }
  }  //   if ( tp == 2 ) ... end-----------------------------------------------



    DP[0] = trialDeform[0] - commitQ[0];
    DP[1] = trialDeform[1] - commitQ[1];
    Beta = sqrt(pow(DP[0],2) + pow(DP[1],2));

    if(Beta > Alpha){
      P_AB = Alpha / Beta;
      trialQ[0] = commitQ[0] *P_AB + trialDeform[0] * (1.0 - P_AB);
      trialQ[1] = commitQ[1] *P_AB + trialDeform[1] * (1.0 - P_AB);
      trialForce[0] = Kxy *(trialDeform[0] - trialQ[0]);
      trialForce[1] = Kxy *(trialDeform[1] - trialQ[1]);
    }else{
      trialQ[0] = commitQ[0];
      trialQ[1] = commitQ[1];
      trialForce[0] = Kxy *(trialDeform[0] - trialQ[0]);
      trialForce[1] = Kxy *(trialDeform[1] - trialQ[1]);
    }


  // 剛性計算
    for (int i=0; i<3; i++) {
      if ((trialDeform[i] - commitDeform[i]) < DBL_EPSILON) {
       trialStiff[i] = initialStiff[i];
      } else {
       trialStiff[i] = initialStiff[i];
      }
    }

    // check print--------------------------------------------
    //opserr << "  mui = " << mui << endln;
    //opserr << "  trialStiff[0] = " << trialStiff[0] << endln;
    //opserr << "  trialStiff[1] = " << trialStiff[1] << endln;
    //opserr << "  trialStiff[2] = " << trialStiff[2] << endln;
    //------------------------------------------------------

  return 0;
}

  const double& EBSBearing::getStrain(int direction)
{
  return trialDeform[direction];
}

  const double& EBSBearing::getStress(int direction)
{
  return trialForce[direction];
}

  const double& EBSBearing::getTangent(int direction)
{
  return trialStiff[direction];
}

  const double& EBSBearing::getInitialTangent(int direction)
{
  return initialStiff[direction];
}


int EBSBearing::commitState()
{
  int errCode = 0;



// commit <- trial

  for (int i=0; i<3; i++) {
    commitDeform[i]  = trialDeform[i];
    commitForce[i]   = trialForce[i];
    commitStiff[i]   = trialStiff[i];
    commitQ[i] = trialQ[i];
  }

  return errCode;
}


int EBSBearing::revertToLastCommit()
{
  int errCode = 0;

  // trial <- commit
  for (int i=0; i<3; i++) {
    trialDeform[i]  = commitDeform[i];
    trialForce[i]   = commitForce[i];
    trialStiff[i]   = commitStiff[i];
    trialQ[i] = commitQ[i];
  }

  return errCode;
}


int EBSBearing::revertToStart()
{
  int errCode=0;

  // trial variables
  basicDisp.Zero();
  basicForce.Zero();

  // reset basic stiffness matrix
  basicStiff = basicStiffInit;

  // データの初期化
  for (int i=0; i<3; i++) {

    trialDeform[i]  = 0.0;
    trialForce[i]   = 0.0;
    trialStiff[i]   = initialStiff[i];
    trialQ[i] = 0.0;

    commitDeform[i] = 0.0;
    commitForce[i] = 0.0;
    commitStiff[i] = initialStiff[i];
    commitQ[i] = 0.0;

  }

  return errCode;
}


int EBSBearing::update()
{
  // get global trial displacements and velocities
  const Vector &dsp1 = theNodes[0]->getTrialDisp();
  const Vector &dsp2 = theNodes[1]->getTrialDisp();
  
  static Vector globalDisp(12), globalDispDot(12);
  for (int i=0; i<6; i++)  {
    globalDisp(i)   = dsp1(i);
    globalDisp(i+6) = dsp2(i);
  }

  static Vector localDispDot(12);

  
  // transform response from the global to the local system
  localDisp    = Tgl*globalDisp;
  
  // transform response from the local to the basic system
  basicDisp    = Tlb*localDisp;


  // calculate shear forces and stiffnesses in basic y- and z-direction
  // get trial shear forces of hysteretic component
  basicForce.Zero();
  basicStiff.Zero();

  this->setTrialStrain(basicDisp);

  basicForce(0) = this->getStress(2);//2014.02.21 EBSBearingの挿入
  basicForce(1) = this->getStress(0);
  basicForce(2) = this->getStress(1);

  basicStiff(0,0) = this->getTangent(2);//2014.02.24 EBSBearingの挿入
  basicStiff(1,1) = this->getTangent(0);
  basicStiff(2,2) = this->getTangent(1);

  return 0;
}


const Matrix& EBSBearing::getTangentStiff()
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


const Matrix& EBSBearing::getInitialStiff()
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


const Matrix& EBSBearing::getMass()
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
    theMatrix(i+3,i+3) = m;
  }
  
  return theMatrix; 
}


void EBSBearing::zeroLoad()
{
  theLoad.Zero();
}

int EBSBearing::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
  opserr <<"EBSBearing::addLoad() - "
	 << "load type unknown for element: "
	 << this->getTag() << endln;
  
  return -1;
}

int EBSBearing::addInertiaLoadToUnbalance(const Vector &accel)
{
  // check for quick return
  if (mass == 0.0)  {
    return 0;
  }
  
  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  
  if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
    opserr << "EBSBearing::addInertiaLoadToUnbalance() - "
	   << "matrix and vector sizes are incompatible\n";
    return -1;
  }
  
  // want to add ( - fact * M R * accel ) to unbalance
  // take advantage of lumped mass matrix
  double m = 0.5*mass;
  for (int i = 0; i < 3; i++)  {
    theLoad(i)   -= m * Raccel1(i);
    theLoad(i+3) -= m * Raccel2(i);
  }
  
  return 0;
}


const Vector& EBSBearing::getResistingForce()
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


const Vector& EBSBearing::getResistingForceIncInertia()
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
      theVector(i+3) += m * accel2(i);
    }
  }
  
  return theVector;
}


int EBSBearing::sendSelf(int commitTag, Channel &sChannel)
{
  return -1;
}


int EBSBearing::recvSelf(int commitTag, Channel &rChannel,
				  FEM_ObjectBroker &theBroker)
{
  return -1;
}


int EBSBearing::displaySelf(Renderer &theViewer,
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


void EBSBearing::Print(OPS_Stream &s, int flag)
{
  
  if (flag == 0)  {
    // print everything
    s << "Element: " << this->getTag(); 
    s << "  type: EBSBearing  iNode: " << connectedExternalNodes(0);
    s << "  jNode: " << connectedExternalNodes(1) << endln;
    
    s << "Input parameters: " << endln;
    s << "  Tp: " << tp << endln;
    s << "  Ar: " << ar << endln;
    s << "  Hr: " << hr << endln;
    s << "  Mu: " << mu << endln;
    s << "  Ec: " << ec << endln;
    s << "  Gr: " << gr << endln;
    s << "  Fz0:" << fz0<< endln;
    s << "  As: " << as << endln;
    
    // determine resisting forces in global system
    s << "  resisting force: " << this->getResistingForce() << endln;
    s << "  Disp-x: " << this->getStrain(0) << endln;
    s << "  Force-x: " << this->getStress(0) << endln;
    s << "  Stiff-x: " << this->getTangent(0) << endln;
    s << "  StiffInit-x: " << this->getInitialTangent(0) << endln;
    s << "  Disp-y: " << this->getStrain(1) << endln;
    s << "  Force-y: " << this->getStress(1) << endln;
    s << "  Stiff-y: " << this->getTangent(1) << endln;
    s << "  StiffInit-y: " << this->getInitialTangent(1) << endln;

  } else if (flag == 1)  {
    // does nothing
  }
}


Response* EBSBearing::setResponse(const char **argv, int argc,
					   OPS_Stream &output)
{
  Response *theResponse = 0;
  
  output.tag("ElementOutput");
  output.attr("eleType","EBSBearing");
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
  
  output.endTag(); // ElementOutput
  
  return theResponse;
}


int EBSBearing::getResponse(int responseID, Information &eleInfo)
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
    
  default:
    return -1;
  }
}


// set up the transformation matrix for orientation
void EBSBearing::setUp()
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
      opserr << "WARNING EBSBearing::setUp() - " 
	     << "element: " << this->getTag() << endln
	     << "ignoring nodes and using specified "
	     << "local x vector to determine orientation\n";
    }
  }
  // check that vectors for orientation are of correct size
  if (oriX.Size() != 3 || oriYp.Size() != 3)  {
    opserr << "EBSBearing::setUp() - "
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
    opserr << "EBSBearing::setUp() - "
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
