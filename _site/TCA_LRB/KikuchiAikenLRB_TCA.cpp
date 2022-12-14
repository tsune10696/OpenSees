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
// $Date: 2015-02-01 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/KikuchiAikenLRB_TCA.cpp,v $

// Written: Ippei Tsunezawa
// Created: February 1, 2015
//
// Kikuchi&Aiken model for lead rubber bearing (for TCA)
//
// Description: This file contains the function to parse the TCL input
// uniaxialMaterial KikuchiAikenLRB_TCA matTag? type? ar? hr? gr? ap? tp? alph? beta? <-T temp? > <-coKQ rk? rq?> <-coMSS rs? rf?>


#include <TclModelBuilder.h>
#include <string.h>
#include <tcl.h>

#include <KikuchiAikenLRB_TCA.h>
#include <Vector.h>
#include <Channel.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>



int
TclCommand_KikuchiAikenLRB_TCA(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  //�Œ���K�v�Ȉ���
  int tag;
  int type = 1;
  double ar = 0.0;
  double hr = 0.0;
  double gr = 0.392e6;
  double ap = 0.0;
  double tp = 8.33e6;
  double alph = 0.588e6;
  double beta = 13.0;

  //�I�v�V���������̃f�t�H���g�l

  double temp = 20.0;
  //double temp = 15.0;
  double rk = 1.0;
  double rq = 1.0;
  double rs = 1.0;
  double rf = 1.0;

  //�쐬����I�u�W�F�N�g
  UniaxialMaterial *theMaterial = 0;


  //���͈����̃G���[�`�F�b�N
  bool ifNoError = true;

  if (argc < 11) { // uniaxialMaterial KikuchiAikenLRB_TCA matTag? type? ar? hr? gr? ap? tp? alph? beta?

    opserr << "WARNING KikuchiAikenLRB_TCA invalid number of arguments\n";
    ifNoError = false;

  } else {

    //argv[2~10]�̃`�F�b�N
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING KikuchiAikenLRB_TCA invalid tag" << endln;
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[3], &type) != TCL_OK || type <= 0 || type >=4 ) { // 2014.06.16 mkiku
      opserr << "WARNING KikuchiAikenLRB_TCA invalid type" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[4], &ar) != TCL_OK || ar <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB_TCA invalid ar" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[5], &hr) != TCL_OK || ar <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB_TCA invalid hr" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[6], &gr) != TCL_OK || gr <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB_TCA invalid gr" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[7], &ap) != TCL_OK || ap <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB_TCA invalid ap" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[8], &tp) != TCL_OK || tp <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB_TCA invalid tp" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[9], &alph) != TCL_OK || alph <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB_TCA invalid alph" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[10], &beta) != TCL_OK || beta <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB_TCA invalid beta" << endln;
      ifNoError = false;
    }


    //argv[11~]�̃`�F�b�N
    for (int i=11; i<=(argc-1); i++) {

      if (strcmp(argv[i],"-T")==0 && (i+1)<=(argc-1)) { // <-T temp?> �̓ǂݍ���

	if (Tcl_GetDouble(interp,argv[i+1], &temp) != TCL_OK) {
	  opserr << "WARNING KikuchiAikenLRB_TCA invalid temp" << endln;
	  ifNoError = false;
	}
	
	i += 1;
	
      } else if (strcmp(argv[i],"-coKQ")==0 && (i+2)<=(argc-1)) { // <-coKQ rk? rq?> �̓ǂݍ���
	
	if (Tcl_GetDouble(interp,argv[i+1], &rk) != TCL_OK || rk < 0.0) {
	  opserr << "WARNING KikuchiAikenLRB_TCA invalid rk" << endln;
	  ifNoError = false;
	}
	
	if (Tcl_GetDouble(interp,argv[i+2], &rq) != TCL_OK || rq < 0.0) {
	  opserr << "WARNING KikuchiAikenLRB_TCA invalid rq" << endln;
	  ifNoError = false;
	} 

	i += 2;
      } else if (strcmp(argv[i],"-coMSS")==0 && (i+2)<=(argc-1)) { // <-coMSS rs? rf?> �̓ǂݍ���
	
	if (Tcl_GetDouble(interp,argv[i+1], &rs) != TCL_OK || rs < 0.0) {
	  opserr << "WARNING KikuchiAikenLRB_TCA invalid rs" << endln;
	  ifNoError = false;
	}
	
	if (Tcl_GetDouble(interp,argv[i+2], &rf) != TCL_OK || rf < 0.0) {
	  opserr << "WARNING KikuchiAikenLRB_TCA invalid rf" << endln;
	  ifNoError = false;
	} 

	i += 2;
	
      } else { // �����ȃI�v�V����
	opserr << "WAINING KikuchiAikenLRB_TCA invalid optional arguments" << endln;
	ifNoError = false;
	break;
      }

    }

  }

  //���͈����̃G���[����
  if (!ifNoError) {
    //���̓f�[�^
    opserr << "Input command: ";
    for (int i=0; i<argc; i++){
      opserr << argv[i] << " ";
    }
    opserr << endln;
    
    //�K�v�f�[�^
    opserr << "Want: uniaxialMaterial KikuchiAikenLRB_TCA matTag? type? ar? hr? gr? ap? tp? alph? beta? <-T temp? > <-coKQ rk? rq?> <-coMSS rs? rf?>" << endln;
//
    return TCL_ERROR;
  }

  //�f�t�H���g�l��ݒ肵�悤�Ƃ��ĊԈႦ���ꍇ�ɑΉ�
  if (rk == 0.0) rk = 1.0;
  if (rq == 0.0) rq = 1.0;
  if (rs == 0.0) rs = 1.0;
  if (rf == 0.0) rf = 1.0;

  // Parsing was successful, allocate the material
  theMaterial = new KikuchiAikenLRB_TCA(tag, type, ar, hr, gr, ap, tp, alph, beta, temp, rk, rq, rs, rf);

  // check print--------------------------------------------
  //  opserr << "   \n";
  //  opserr << "TclCommand_KikuchiAikenLRB_TCA\n";
  //  opserr << "  tag  = " << tag << endln;
  //  opserr << "  type = " << type << endln;
  //  opserr << "  ar   = " << ar << endln;
  //  opserr << "  hr   = " << hr << endln;
  //  opserr << "  gr   = " << gr << endln;
  //  opserr << "  ap   = " << ap << endln;
  //  opserr << "  tp   = " << tp << endln;
  //  opserr << "  alph = " << alph << endln;
  //  opserr << "  beta = " << beta << endln;
  //  opserr << "  temp = " << temp << endln;
  //  opserr << "  rk   = " << rk << endln;
  //  opserr << "  rq   = " << rq << endln;
  //  opserr << "  rs   = " << rs << endln;
  //  opserr << "  rf   = " << rf << endln;
  // -------------------------------------------------------

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial " << argv[1] << endln;
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (OPS_addUniaxialMaterial(theMaterial) == false) {
    opserr << "WARNING could not add uniaxialMaterial to the modelbuilder\n";
    opserr << *theMaterial << endln;
    delete theMaterial; // invoke the material objects destructor, otherwise mem leak
    return TCL_ERROR;
  } else {
    return TCL_OK;
  }

}




KikuchiAikenLRB_TCA::KikuchiAikenLRB_TCA(int tag, int type, double ar, double hr, double gr, double ap, double tp, 
		  double alph, double beta, double temp, double rk, double rq, double rs, double rf)
  :UniaxialMaterial(tag,MAT_TAG_KikuchiAikenLRB_TCA),Type(type),Ar(ar),Hr(hr),Gr(gr),Ap(ap),Tp(tp),
   Alph(alph),Beta(beta),Temp(temp),Rk(rk),Rq(rq),Rs(rs),Rf(rf)
{

  //�S����ɉ������e��p�����[�^�̐ݒ�E�֐��̑I��
  switch (Type) {

  case 1: // Oiles LRB 500mm�a�ی^ S2=4, �ʈ�1MPa (idac OILES500F����)
    trgStrain = 0.05;
    lmtStrain = 4.10;
    calcN   = KikuchiAikenLRB_TCA::calcNType1;
    calcP   = KikuchiAikenLRB_TCA::calcPType1;
    calcA   = KikuchiAikenLRB_TCA::calcAType1;
    calcB   = KikuchiAikenLRB_TCA::calcBType1;
    calcC   = KikuchiAikenLRB_TCA::calcCType1;
    calcCQd = KikuchiAikenLRB_TCA::calcCQdType1;
    calcCKd = KikuchiAikenLRB_TCA::calcCKdType1;
    calcCHeq= KikuchiAikenLRB_TCA::calcCHeqType1;
    break;

  case 2: // Oiles LRB 250mm�p�^ S2=5, �ʈ�1MPa (idac OILES250S����)
    trgStrain = 0.05;
    lmtStrain = 4.10;
    calcN   = KikuchiAikenLRB_TCA::calcNType2;
    calcP   = KikuchiAikenLRB_TCA::calcPType2;
    calcA   = KikuchiAikenLRB_TCA::calcAType2;
    calcB   = KikuchiAikenLRB_TCA::calcBType2;
    calcC   = KikuchiAikenLRB_TCA::calcCType2;
    calcCQd = KikuchiAikenLRB_TCA::calcCQdType2;
    calcCKd = KikuchiAikenLRB_TCA::calcCKdType2;
    calcCHeq= KikuchiAikenLRB_TCA::calcCHeqType2;
    break;

  case 3: // Standard 400% (idac STARNDARD-500-2����), 2014.06.16 mkiku
    trgStrain = 0.05;
    lmtStrain = 4.10;
    calcN   = KikuchiAikenLRB_TCA::calcNType3;
    calcP   = KikuchiAikenLRB_TCA::calcPType3;
    calcA   = KikuchiAikenLRB_TCA::calcAType3;
    calcB   = KikuchiAikenLRB_TCA::calcBType3;
    calcC   = KikuchiAikenLRB_TCA::calcCType3;
    calcCQd = KikuchiAikenLRB_TCA::calcCQdType3;
    calcCKd = KikuchiAikenLRB_TCA::calcCKdType3;
    calcCHeq= KikuchiAikenLRB_TCA::calcCHeqType3;
    break;
  }

  // �����v�Z
//-------------------------------------------
  Alph_T = 0.4 + 0.25*(Temp / 327.5);
  Tp_TCA = 15.0e6*(1.0 - pow(Temp / 327.5, Alph_T));
  //Tp_ini = Tp_TCA ;
  Tp_ini = 8.33e6 ;
  //ce = Tp_TCA / Tp_ini;


  //qd100 = Tp_TCA*Ap * Rq;                   // �ؕЉ׏d�̊�l
  qd100 = Tp_ini*Ap * Rq;                   // �ؕЉ׏d�̊�l
  kd100 = (Gr*Ar/Hr + Alph*Ap/Hr) * Rk; // �~���㍄���̊�l
  ku100 = Beta * kd100; // ���������̊�l
//-------------------------------------------

  // ���������̌v�Z
  qd = (this->calcCQd)(trgStrain)*qd100;
  kd = (this->calcCKd)(trgStrain)*kd100;
  ku = (this->calcCKd)(trgStrain)*ku100;
  initialStiff = compKeq(fabs(trgStrain*Hr),qd,kd);

  // check print--------------------------------------------
  //  opserr << "   \n";
  //  opserr << "KikuchiAikenLRB_TCA::KikuchiAikenLRB_TCA()\n";
  //  opserr << "  qd100 = " << qd100 << endln;
  //  opserr << "  kd100 = " << kd100 << endln;
  //  opserr << "  ku100 = " << ku100 << endln;
  //  opserr << "  initialStiff = " << initialStiff << endln;
  // -------------------------------------------------------

  //�����Ȑ��L���p�ϐ�
  numIdx = 500; //����500��
  revXBgn   = new double [numIdx];
  revQ2Bgn  = new double [numIdx];
  revXEnd   = new double [numIdx];
  revQ2End  = new double [numIdx];
  revB      = new double [numIdx];
  revAlpha  = new double [numIdx];

  trialDeform  = 0.0;
  trialForce   = 0.0;
  trialStiff   = initialStiff;
  trialStrain  = 0.0;
  trialIfElastic = true;
  trialQ1 = 0.0;
  trialQ2 = 0.0;
  trialMaxStrain = 0.0;
  trialDDeform = 0.0;
  trialDDeformLastSign = 0;
  trialIdxRev=0;

  commitDeform  = 0.0;
  commitForce   = 0.0;
  commitStiff   = initialStiff;
  commitStrain  = 0.0;
  commitIfElastic = true;
  commitQ1 = 0.0;
  commitQ2 = 0.0;
  commitMaxStrain = 0.0;
  commitDDeform = 0.0;
  commitDDeformLastSign = 0;
  commitIdxRev=0;

  revB[0] = 0.0;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenLRB_TCA::KikuchiAikenLRB_TCA\n";
  // -------------------------------------------------------

}

KikuchiAikenLRB_TCA::KikuchiAikenLRB_TCA()
  :UniaxialMaterial(0,MAT_TAG_KikuchiAikenLRB_TCA)
{

  //�����Ȑ��L���p�ϐ�?

  trialDeform  = 0.0;
  trialForce   = 0.0;
  trialStiff   = 0.0;
  commitDeform = 0.0;
  commitForce  = 0.0;
  commitStiff  = 0.0;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenLRB_TCA::KikuchiAikenLRB_TCA()\n";
  // -------------------------------------------------------

}

KikuchiAikenLRB_TCA::~KikuchiAikenLRB_TCA()
{

  if (revXBgn != 0)
    delete [] revXBgn;

  if (revQ2Bgn != 0)
    delete [] revQ2Bgn;

  if (revXEnd != 0)
    delete [] revXEnd;

  if (revQ2End != 0)
    delete [] revQ2End;

  if (revB != 0)
    delete [] revB;

  if (revAlpha != 0)
    delete [] revAlpha;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenLRB_TCA::~KikuchiAikenLRB_TCA\n";
  // -------------------------------------------------------

}

int 
KikuchiAikenLRB_TCA::setTrialStrain(double strain, double strainRate)
{

  // ���𑥂����̎菇�Ōv�Z����
  // (����)       (�o��)
  //  �ψ� ---->   �׏d
  //       ----> �΂˒萔


  //(����: �ψ�->�Ђ���)
  trialDeform = strain;
  trialDeform = trialDeform * (Rs/Rf); //for MSS model
  trialStrain = trialDeform/Hr;


  //�ψʑ���
  trialDDeform = trialDeform - commitDeform;


  // check print--------------------------------------------
  //  opserr << "   \n";
  //  opserr << "KikuchiAikenLRB_TCA::setTrialStrain\n";
  //  opserr << "  trialDeform  = " << trialDeform << endln;
  //  opserr << "  trialDDeform = " << trialDDeform << endln;
  //  opserr << "  trialStrain  = " << trialStrain << endln;
  // -------------------------------------------------------


  //�ψʑ����̂Ȃ��v�Z�X�e�b�v
  if ( fabs(trialDDeform) < DBL_EPSILON ) { // if ( trialDDeform == 0.0 )

//------------------------------------------------------------
    Alph_T = 0.4 + 0.25*(Temp / 327.5);
    Tp_TCA = 15.0e6*(1.0 - pow(Temp / 327.5, Alph_T));

    ce = Tp_TCA / Tp_ini;

//------------------------------------------------------------

    trialForce   = commitQ1    + commitQ2    * ce;
    trialStiff   = commitQ1Stf + commitQ2Stf * ce;

    return 0;
  }


  //�ψʑ����̕��� ���X�e�b�v�p
  if (trialDDeform > 0) {
    trialDDeformLastSign = +1;
  } else if (trialDDeform < 0) {
    trialDDeformLastSign = -1;
  } else {
    trialDDeformLastSign = commitDDeformLastSign;
  }

  //���]�|�C���g�̎Q��
  trialIdxRev = commitIdxRev;


  //�K�p���E�̃`�F�b�N
  if ( fabs(trialStrain) > lmtStrain ) {
    opserr << "uniaxialMaterial KikuchiAikenLRB_TCA: \n";
    opserr << "   Response value exceeded limited strain.\n";
    //return -1;
  }

  //�e�����E�̃`�F�b�N
  if ( fabs(trialStrain) > trgStrain ) {
      trialIfElastic = false;
  }

  //�ő�Ђ��݂̃`�F�b�N�E�X�V
  if ( fabs(trialStrain) > commitMaxStrain ) {
    trialMaxStrain = fabs(trialStrain);
  }

  //�p�����[�^���Z�o
  //�e���͈͂ł͌��݂̕ψʂɉ�����Xm,Fm���X�V
  //�e�Y���͈͂ł͍ő�ψʍX�V���Ɋe�p�����[�^���X�V(�X�P���g���J�[�u��)

  if ( trialIfElastic || fabs(trialStrain) == trialMaxStrain ) {

    //�e���͈͂ł�CQd,CKd,Cheq���̎Z��ɒe�����E�Ђ��݁A�e�����E�ό`��p����
    tmpStrain = (fabs(trialStrain)>trgStrain) ? fabs(trialStrain) : trgStrain ; //max(fabs(trialStrain),trgStrain)
    tmpDeform = tmpStrain * Hr;

    qd = (this->calcCQd)(tmpStrain)*qd100;
    kd = (this->calcCKd)(tmpStrain)*kd100;
    ku = (this->calcCKd)(tmpStrain)*ku100;

    keq = compKeq(tmpDeform,qd,kd); // ��������
    heq = (this->calcCHeq)(tmpStrain) * compHeq(tmpDeform,qd,kd,ku); // �����S�������萔
    u = qd / (keq*tmpDeform); // �׏d�ؕД�

    xm  = fabs(trialDeform);
    fm  = keq*xm;

    n = (this->calcN)(fabs(trialStrain));
    p = (this->calcP)(fabs(trialStrain));
    c = (this->calcC)(fabs(trialStrain));
    a = (this->calcA)(fabs(trialStrain),heq,u);

  }  // if ( trialIfElastic || fabs(trialStrain) == trialMaxStrain )  ... end


  // ����ψ�x
  x = (xm>0.0) ? trialDeform/xm : 0.0;

  // ���]�̏���
  if (!trialIfElastic) {

    //���ׂ܂��͔��]
    if (trialDDeform*commitDDeformLastSign < 0) {


      if ( trialIdxRev == 0 ) { //�X�P���g���J�[�u����̏��ׂ��J�n����ꍇ

	trialIdxRev = 1;

	revXBgn[1]  = commitDeform/xm;
	revQ2Bgn[1] = commitQ2;
	
	//b
	double tqd = (this->calcCQd)(fabs(commitStrain))*qd100;
	double tkd = (this->calcCKd)(fabs(commitStrain))*kd100;
	double tku = (this->calcCKd)(fabs(commitStrain))*ku100;
	double tkeq = compKeq(fabs(commitDeform),tqd,tkd);
	double theq = (this->calcCHeq)(fabs(commitStrain)) * compHeq(fabs(commitDeform),tqd,tkd,tku);
	double tu = tqd/(tkeq*fabs(commitDeform));
	b = (this->calcB)(fabs(commitStrain),a,c,theq,tu);

	revB[1] = b;
	revAlpha[1] = 1.0;


      } else { //���]�̏ꍇ

	trialIdxRev++;


	if (trialIdxRev >= numIdx ) { //�����Ȑ��L���p�ϐ��̗̈�𓮓I�Ɋm�ۂ���
	  int newIdx;
	  double *newXBgn;
	  double *newQ2Bgn;
	  double *newXEnd;
	  double *newQ2End;
	  double *newB;
	  double *newAlpha;
	  
	  newIdx = numIdx + 500; //500���ǉ�
	  newXBgn  = new double [newIdx];
	  newQ2Bgn = new double [newIdx];
	  newXEnd  = new double [newIdx];
	  newQ2End = new double [newIdx];
	  newB     = new double [newIdx];
	  newAlpha = new double [newIdx];
	  
	  for(int i = 0; i<numIdx; i++){
	    newXBgn[i]  = revXBgn[i];
	    newQ2Bgn[i] = revQ2Bgn[i];
	    newXEnd[i]  = revXEnd[i];
	    newQ2End[i] = revQ2End[i];
	    newB[i]     = revB[i];
	    newAlpha[i] = revAlpha[i];
	  }
	  
	  numIdx = newIdx;

	  delete [] revXBgn;
	  delete [] revQ2Bgn;
	  delete [] revXEnd;
	  delete [] revQ2End;
	  delete [] revB;
	  delete [] revAlpha;

	  revXBgn  = newXBgn;
	  revQ2Bgn = newQ2Bgn;
	  revXEnd  = newXEnd;
	  revQ2End = newQ2End;
	  revB     = newB;
	  revAlpha = newAlpha;
	}


	revXEnd[trialIdxRev]  = revXBgn[trialIdxRev-1];
	revQ2End[trialIdxRev] = revQ2Bgn[trialIdxRev-1];
	revXBgn[trialIdxRev]  = commitDeform/xm;
	revQ2Bgn[trialIdxRev] = commitQ2;
	
	//b
	if ( revB[trialIdxRev-1] == 0 ) { // b=0�ł̏��ׂ܂��͔��]�ȍ~
	  b = 0;
	} else if (revXEnd[trialIdxRev]*revXBgn[trialIdxRev] > 0) { //�������̈�Ŕ��]
	  b = 0;
	} else { //�ٕ����̈�Ŕ��]  b�������_�ł̂Ђ��݂ōĕ]��
	  double tqd = (this->calcCQd)(fabs(commitStrain))*qd100;
	  double tkd = (this->calcCKd)(fabs(commitStrain))*kd100;
	  double tku = (this->calcCKd)(fabs(commitStrain))*ku100;
	  double tkeq = compKeq(fabs(commitDeform),tqd,tkd);
	  double theq = (this->calcCHeq)(fabs(commitStrain)) * compHeq(fabs(commitDeform),tqd,tkd,tku);
	  double tu = tqd/(tkeq*fabs(commitDeform)); 
	  b = (this->calcB)(fabs(commitStrain),a,c,theq,tu);
	}

	
	//alpha
	if (trialDDeform > 0) {
	  alpha = (this->compAlpha)(a,revB[trialIdxRev-1],b,c,revXEnd[trialIdxRev],revXBgn[trialIdxRev],revAlpha[trialIdxRev-1]);
	} else {
	  alpha = (this->compAlpha)(a,revB[trialIdxRev-1],b,c,-revXEnd[trialIdxRev],-revXBgn[trialIdxRev],revAlpha[trialIdxRev-1]);
	}

	revB[trialIdxRev] = b;
	revAlpha[trialIdxRev] = alpha;

      }

    }


    //���]�|�C���g�̉��� 
    if  (fabs(trialStrain) == trialMaxStrain) { //�X�P���g���J�[�u��,���ׂĉ���
      trialIdxRev = 0;
    } else if (trialIdxRev >= 2 ) { //�J�n�_�Ǝw���_�̊Ԃɕψʓ_���Ȃ��Ƃ�����
      while (trialIdxRev >= 2 && (x-revXBgn[trialIdxRev])*(x-revXEnd[trialIdxRev]) > 0) {
	trialIdxRev--;
      }
    }

  }  //   if (!trialIfElastic)  ... end

  // ���͂̎Z�o

  // Q1 (����)
  if (trialStrain > 0) {
    trialQ1 = (this->compQ1)(u,n,p,fm,x);
    //    q1Stf   = (this->compQ1Derivertive)(u,n,p,keq,x);
    trialQ1Stf   = (this->compQ1Derivertive)(u,n,p,keq,x);
  } else {
    trialQ1 = (this->compQ1)(u,n,p,-fm,-x);
    //    q1Stf   = (this->compQ1Derivertive)(u,n,p,keq,-x);
    trialQ1Stf   = (this->compQ1Derivertive)(u,n,p,keq,-x);
  }

  //Q2(����)
  if (trialIdxRev == 0) {// �X�P���g���J�[�u��

    if (trialDDeform > 0) {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,-fm,-x);
      //q2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,x);
      trialQ2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,x);
    } else {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,fm,x);
      //q2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,-x);
      trialQ2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,-x);
    }

  } else if (trialIdxRev == 1) { //�X�P���g���J�[�u����̏���

    if (trialDDeform > 0) {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,fm,x);
      //      q2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,x);
      trialQ2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,x);
    } else {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,-fm,-x);
      //      q2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,-x);
      trialQ2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,-x);
    }

  } else { //���]��

    if (trialDDeform > 0) {
      trialQ2 = (this->compQ2Masing)(u,a,revB[trialIdxRev],c,fm,x,revXBgn[trialIdxRev],revQ2Bgn[trialIdxRev],revAlpha[trialIdxRev]);
      //q2Stf   = (this->compQ2MasingDerivertive)(u,a,revB[trialIdxRev],c,keq,x,revXBgn[trialIdxRev],revAlpha[trialIdxRev]);
      trialQ2Stf   = (this->compQ2MasingDerivertive)(u,a,revB[trialIdxRev],c,keq,x,revXBgn[trialIdxRev],revAlpha[trialIdxRev]);
    } else {
      trialQ2 = (this->compQ2Masing)(u,a,revB[trialIdxRev],c,-fm,-x,-revXBgn[trialIdxRev],revQ2Bgn[trialIdxRev],revAlpha[trialIdxRev]);
      //q2Stf   = (this->compQ2MasingDerivertive)(u,a,revB[trialIdxRev],c,keq,-x,-revXBgn[trialIdxRev],revAlpha[trialIdxRev]);
      trialQ2Stf   = (this->compQ2MasingDerivertive)(u,a,revB[trialIdxRev],c,keq,-x,-revXBgn[trialIdxRev],revAlpha[trialIdxRev]);
    }

  }

//------------------------------------------------------------
  Alph_T = 0.4 + 0.25*(Temp / 327.5);
  Tp_TCA = 15.0e6*(1.0 - pow(Temp / 327.5, Alph_T));
  ce = Tp_TCA / Tp_ini;
//------------------------------------------------------------


  //�o�́F�׏d
  trialForce  = trialQ1 + trialQ2 * ce;

  //�o�́F����
  if (trialIfElastic) {
    trialStiff = initialStiff;
  } else {
    trialStiff = trialQ1Stf + trialQ2Stf * ce;
  }

  trialForce = trialForce * Rf; //for MSS model
  trialStiff = trialStiff * Rs; //for MSS model


  //---�ǉ�---------------------------
  trialQ1 = trialQ1 * Rf; //for MSS model
  trialQ1Stf = trialQ1Stf * Rs; //for MSS model

  trialQ2 = trialQ2 * Rf; //for MSS model
  trialQ2Stf = trialQ2Stf * Rs; //for MSS model
  //----------------------------------

  return 0;

}

int 
KikuchiAikenLRB_TCA::setTemp(double target_temp){
  Temp = target_temp;
  return 0;
}

double
KikuchiAikenLRB_TCA::getTemp(void){
  return Temp;

}


double 
KikuchiAikenLRB_TCA::getStress(void)
{
  return trialForce;
}


//--�ǉ�-----------------
double 
KikuchiAikenLRB_TCA::getQ1(void)
{
  return trialQ1;
}

double 
KikuchiAikenLRB_TCA::getQ2(void)
{
  return trialQ2 * ce;
}

double 
KikuchiAikenLRB_TCA::getHr(void)
{
  return Hr;
}

double 
KikuchiAikenLRB_TCA::getAr(void)
{
  return Ar;
}

double 
KikuchiAikenLRB_TCA::getAp(void)
{
  return Ap;
}

double 
KikuchiAikenLRB_TCA::getQd(void)
{
  return Tp_TCA * Ap;
}

//--�ǉ�-----------------


double 
KikuchiAikenLRB_TCA::getTangent(void)
{

  return trialStiff;
}

double 
KikuchiAikenLRB_TCA::getInitialTangent(void)
{
  return initialStiff;
}

double 
KikuchiAikenLRB_TCA::getStrain(void)
{
  return trialDeform;
}

int 
KikuchiAikenLRB_TCA::commitState(void)
{
  commitDeform  = trialDeform;
  commitForce   = trialForce;
  commitStiff   = trialStiff;
  commitStrain  = trialStrain;
  commitIfElastic = trialIfElastic;
  commitQ1 = trialQ1;
  commitQ2 = trialQ2;
  commitMaxStrain = trialMaxStrain;
  commitDDeform   = trialDDeform;
  commitDDeformLastSign = trialDDeformLastSign;
  commitIdxRev = trialIdxRev;
  commitQ1Stf = trialQ1Stf;
  commitQ2Stf = trialQ2Stf;

  //opserr << "ce=" << ce << " Tp_ini=" << Tp_ini << " Tp_TCA=" << Tp_TCA << endln;

  //opserr << Tp_TCA  << endln;

  
  

  // check print--------------------------------------------
  // opserr << "KikuchiAikenLRB_TCA::commitState\n";
  // -------------------------------------------------------
  return 0;
}

int 
KikuchiAikenLRB_TCA::revertToLastCommit(void)
{

  trialDeform  = commitDeform;
  trialForce   = commitForce;
  trialStiff   = commitStiff;
  trialStrain  = commitStrain;
  trialIfElastic = commitIfElastic;
  trialQ1 = commitQ1;
  trialQ2 = commitQ2;
  trialMaxStrain = commitMaxStrain;
  trialDDeform = commitDDeform;
  trialDDeformLastSign = commitDDeformLastSign;
  trialIdxRev = commitIdxRev;
  trialQ1Stf = commitQ1Stf;
  trialQ2Stf = commitQ2Stf;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenLRB_TCA::revertToLastCommit\n";
  // -------------------------------------------------------
  return 0;
}

int 
KikuchiAikenLRB_TCA::revertToStart(void)
{

  trialDeform  = 0.0;
  trialForce   = 0.0;
  trialStiff   = initialStiff;
  trialStrain  = 0.0;
  trialIfElastic = true;
  trialQ1 = 0.0;
  trialQ2 = 0.0;
  trialMaxStrain = 0.0;
  trialDDeform = 0.0;
  trialDDeformLastSign = 0;
  trialIdxRev=0;

  commitDeform  = 0.0;
  commitForce   = 0.0;
  commitStiff   = initialStiff;
  commitStrain  = 0.0;
  commitIfElastic = true;
  commitQ1 = 0.0;
  commitQ2 = 0.0;
  commitMaxStrain = 0.0;
  commitDDeform = 0.0;
  commitDDeformLastSign = 0;
  commitIdxRev=0;



  //----------------------------

    trialQ1Stf  = initialStiff  / 2.0;
    commitQ1Stf = initialStiff  / 2.0;
    trialQ2Stf  = initialStiff  / 2.0;
    commitQ2Stf = initialStiff  / 2.0;
  //----------------------------


  revB[0] = 0.0;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenLRB_TCA::revertToStart\n";
  // -------------------------------------------------------
  return 0;
}

UniaxialMaterial *
KikuchiAikenLRB_TCA::getCopy(void)  //�����Ȑ��L���p�ϐ��̎󂯓n��?
{

  KikuchiAikenLRB_TCA *theCopy = new KikuchiAikenLRB_TCA(this->getTag(), Type, Ar, Hr, Gr, Ap, Tp, 
			      Alph, Beta, Temp, Rk, Rq, Rs, Rf );

  //  theCopy->qd100 = qd100;
  //  theCopy->kd100 = kd100;
  //  theCopy->ku100 = ku100;
  //  theCopy->qd = qd;
  //  theCopy->kd = kd;
  //  theCopy->ku = ku;

  //  theCopy->trgStrain = trgStrain;
  //  theCopy->lmtStrain = lmtStrain;
  //  theCopy->initialStiff = initialStiff;

  // Copy temporary variables
  theCopy->tmpStrain = tmpStrain;
  theCopy->keq = keq;
  theCopy->heq = heq;
  theCopy->u = u;
  theCopy->n = n;
  theCopy->a = a;
  theCopy->b = b;
  theCopy->c = c;
  theCopy->xm = xm;
  theCopy->fm = fm;
  theCopy->x = x;
  theCopy->alpha = alpha;

  // Copy trial variables
  theCopy->trialDeform = trialDeform;
  theCopy->trialForce = trialForce;
  theCopy->trialStiff = trialStiff;
  theCopy->trialStrain = trialStrain;
  theCopy->trialIfElastic = trialIfElastic;
  theCopy->trialQ1 = trialQ1;
  theCopy->trialQ2 = trialQ2;
  theCopy->trialMaxStrain = trialMaxStrain;
  theCopy->trialDDeform = trialDDeform;
  theCopy->trialDDeformLastSign = trialDDeformLastSign;
  theCopy->trialIdxRev = trialIdxRev;

  // Copy commit variables
  theCopy->commitDeform = commitDeform;
  theCopy->commitForce = commitForce;
  theCopy->commitStiff = commitStiff;
  theCopy->commitStrain = commitStrain;
  theCopy->commitIfElastic = commitIfElastic;
  theCopy->commitQ1 = commitQ1;
  theCopy->commitQ2 = commitQ2;
  theCopy->commitMaxStrain = commitMaxStrain;
  theCopy->commitDDeform = commitDDeform;
  theCopy->commitDDeformLastSign = commitDDeformLastSign;
  theCopy->commitIdxRev = commitIdxRev;

  // check print--------------------------------------------
  // opserr << "KikuchiAikenLRB_TCA::getCopy\n";
  // -------------------------------------------------------

  return theCopy;
}

int 
KikuchiAikenLRB_TCA::sendSelf(int cTag, Channel &theChannel) //�����Ȑ��L���p�ϐ��̎󂯓n��?
{

//   int res = 0;

//   static Vector data(51);

//   data(0)  = this->getTag();
//   data(1)  = Type; //int -> double
//   data(2)  = Ar;
//   data(3)  = Hr;
//   data(4)  = Gr;
//   data(5)  = Ap;
//   data(6)  = Tp;
//   data(7)  = Alph;
//   data(8)  = Beta;
//   data(9)  = Temp;
//   data(10) = Rk;
//   data(11) = Rq;
//   data(12) = Rs;
//   data(13) = Rf;
//   data(14  = trgStrain;
//   data(15) = lmtStrain;
//   data(16) = initialStiff;
//   data(17) = tmpstrain;
//   data(18) = keq;
//   data(19) = heq;
//   data(20) = u;
//   data(21) = n;
//   data(22) = p;
//   data(23) = a;
//   data(24) = b;
//   data(25) = c;
//   data(26) = xm;
//   data(27) = fm;
//   data(28) = x;
//   data(29) = alpha;
//   data(30) = trialDeform;
//   data(31) = trialForce;
//   data(32) = trialStiff;
//   data(33) = trialStrain;
//   data(34) = (trialIfElastic) ? 1 : 0; //bool -> double
//   data(35) = trialQ1;
//   data(36) = trialQ2;
//   data(37) = trialMaxStrain;
//   data(38) = trialDDeform;
//   data(39) = trialDDeformLastSign; //int -> double
//   data(40) = trialIdxRev; //int -> double
//   data(41) = commitDeform;
//   data(42) = commitForce;
//   data(43) = commitStiff;
//   data(44) = commitStrain;
//   data(45) = (commitIfElastic) ? 1 : 0; //bool -> double
//   data(46) = commitQ1;
//   data(47) = commitQ2;
//   data(48) = commitMaxStrain;
//   data(49) = commitDDeform;
//   data(50) = commitDDeformLastSign; //int -> double
//   data(51) = commitIdxRev; //int -> double

//   res = theChannel.sendVector(this->getDbTag(), cTag, data);
//   if (res < 0) 
//     opserr << "KikuchiAikenLRB_TCA::sendSelf() - failed to send data\n";

//   return res;
  
  return -1;

}

int 
KikuchiAikenLRB_TCA::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker) //�����Ȑ��L���p�ϐ��̎󂯓n��?
{
//   int res = 0;

//   static Vector data(51);
//   res = theChannel.recvVector(this->getDbTag(), cTag, data);

//   if (res < 0) {
//       opserr << "KikuchiAikenLRB_TCA::recvSelf() - failed to receive data\n";
//       this->setTag(0);      
//   }
//   else {
//     this->setTag((int)data(0));
//     Type = (int)data(1); //double -> int
//     Ar = data(2);
//     Hr = data(3);
//     Gr = data(4);
//     Ap = data(5);
//     Tp = data(6);
//     Alph = data(7);
//     Beta = data(8);
//     Temp = data(9);
//     Rk = data(10);
//     Rq = data(11);
//     Rs = data(12);
//     Rf = data(13);
//     trgStrain = data(14);
//     lmtStrain = data(15);
//     initialStiff = data(16);
//     tmpstrain = data(17);
//     keq = data(18);
//     heq = data(19);
//     u = data(20);
//     n = data(21);
//     p = data(22);
//     a = data(23);
//     b = data(24);
//     c = data(25);
//     xm = data(26);
//     fm = data(27);
//     x = data(28);
//     alpha = data(29);
//     trialDeform = data(30);
//     trialForce = data(31);
//     trialStiff = data(32);
//     trialStrain = data(33);
//     trialIfElastic = ( (int)data(34) == 1) ? true : false; //double -> bool
//     trialQ1 = data(35);
//     trialQ2 = data(36);
//     trialMaxStrain = data(37);
//     trialDDeform = data(38);
//     trialDDeformLastSign = (int)data(39); //int -> double
//     trialIdxRev = (int)data(40); //int -> double
//     commitDeform = data(41);
//     commitForce = data(42);
//     commitStiff = data(43);
//     commitStrain = data(44);
//     commitIfElastic = ( (int)data(45) == 1) ? true : false; //double -> bool
//     commitQ1 = data(46);
//     commitQ2 = data(47);
//     commitMaxStrain = data(48);
//     commitDDeform = data(49);
//     commitDDeformLastSign = (int)data(50); //int -> double
//     commitIdxRev = (int)data(51); //int -> double
//   }

//   return res;

  return -1;
}

void 
KikuchiAikenLRB_TCA::Print(OPS_Stream &s, int flag)
{
  s << "KikuchiAikenLRB_TCA : " << this->getTag() << endln;
  s << "  Type: " << Type << endln;
  s << "  Ar: " << Ar << endln;
  s << "  Hr: " << Hr << endln;
  s << "  Gr: " << Gr << endln;
  s << "  Ap: " << Ap << endln;
  s << "  Tp: " << Tp << endln;
  s << "  Alph: " << Alph << endln;
  s << "  Beta: " << Beta << endln;
  s << "  Temp: " << Temp << endln;
  s << "  Rk: " << Rk << endln;
  s << "  Rq: " << Rq << endln;
  s << "  Rs: " << Rs << endln;
  s << "  Rf: " << Rf << endln;

  return;
}



//------------------------------------------------------------
//�e�n���f���v�Z�p�̊֐�
//------------------------------------------------------------

double
KikuchiAikenLRB_TCA::compQ1(double u, double n, double p, double fm, double x)
{
//  return 0.5*(1-u)*fm*(x+pow(x,n));
//  new F1 model, u�̕����ɂ��ẮA�֐��̌Ăяo�����ɍl��
  return (1-u)*fm*( (1-p)*x + p*pow(x,n) );
}

double
KikuchiAikenLRB_TCA::compQ2Unload(double u, double a, double b, double c, double fm, double x)
{
  return u*fm*(1-2*exp(-a*(1+x))+b*(1+x)*exp(-c*(1+x)));
}

double
KikuchiAikenLRB_TCA::compQ2Masing(double u, double a, double b, double c, double fm, double x1, double x2, double q2i, double alpha)
{
  return q2i + alpha*u*fm*(2-2*exp(-a*(x1-x2))+b*(x1-x2)*exp(-c*(x1-x2)));
}


double
KikuchiAikenLRB_TCA::compAlpha(double a, double b1, double b2, double c, double x1, double x2, double alpha0)
{
  return alpha0*(2-2*exp(-a*(x1-x2))+b1*(x1-x2)*exp(-c*(x1-x2)))/(2-2*exp(-a*(x1-x2))+b2*(x1-x2)*exp(-c*(x1-x2)));
}


double
KikuchiAikenLRB_TCA::compQ1Derivertive(double u, double n, double p, double keq, double x)
{
//  return 0.5*(1-u)*keq*(1+n*pow(x,n-1));
//  new F1 model, u�̕����ɂ��ẮA�֐��̌Ăяo�����ɍl��
    return (1-u)*keq*( (1-p) + p*n*pow(x,n-1) );
}

double
KikuchiAikenLRB_TCA::compQ2UnloadDerivertive(double u, double a, double b, double c, double keq, double x)
{
  return u*keq*(2*a*exp(-a*(1+x))+b*exp(-c*(1+x))-b*c*(1+x)*exp(-c*(1+x)));
}


double
KikuchiAikenLRB_TCA::compQ2MasingDerivertive(double u, double a, double b, double c, double keq, double x1, double x2,double alpha)
{
  return alpha*u*keq*(2*a*exp(-a*(x1-x2))+b*exp(-c*(x1-x2))-b*c*(x1-x2)*exp(-c*(x1-x2)));
}

//�p�����[�^a�Z�o�p�̓񕪖@
double
KikuchiAikenLRB_TCA::compABisection(double heq, double u, double min, double max, double tol, double lim)
{
  double aMin, aMax, aTmp ;
  double RHS,LHS ;

  RHS = (2*u-M_PI*heq)/(2*u);

  aMin = min;
  aMax = max;

  while (true) {
    aTmp = (aMin+aMax)/2;
    LHS = (1-exp(-2*aTmp))/(aTmp);
    if (fabs((LHS-RHS)/RHS) < tol) {
      break;
    } else if (LHS < RHS) {
      aMax = aTmp;
    } else {
      aMin = aTmp;
    }
  }

  return (aTmp < lim) ? aTmp : lim ; //min(aTmp,lim)

}

//���������̌v�Z
double
KikuchiAikenLRB_TCA::compKeq(double xm, double qd, double kd)
{
  return (qd + kd*xm)/xm;
}

//�����S�������萔�̌v�Z
double
KikuchiAikenLRB_TCA::compHeq(double xm, double qd, double kd, double ku)
{
  double dy,qy,fm,heq;

    dy  = qd/(ku-kd);
    qy  = ku*dy;
    fm  = qd + kd*xm;
    heq = 2.0*qd/M_PI/fm*(1.0-qy/ku/xm);

  return heq ;

}



//------------------------------------------------------------
//�S����ɉ������e��֐��̐ݒ�
//------------------------------------------------------------

//--------------------------
// type1 : OILES LRB500R
//--------------------------

double KikuchiAikenLRB_TCA::calcNType1(double gm)
{return 8.0;}

double KikuchiAikenLRB_TCA::calcPType1(double gm)
{
  if      (gm<2.0)  {return 0.0;}
  else              {return -0.30226 + 0.15113*gm;}
}

double KikuchiAikenLRB_TCA::calcCType1(double gm)
{return 6.0;}

double KikuchiAikenLRB_TCA::calcCQdType1(double gm)
{
  if      (gm<0.1)  {return 2.036*pow(gm,0.410);}
  else if (gm<0.5)  {return 1.106*pow(gm,0.145);}
  else              {return 1.0;}
}

double KikuchiAikenLRB_TCA::calcCKdType1(double gm)
{
  if      (gm<0.25) {return 0.779*pow(gm,-0.43 );}
  else if (gm<1.0)  {return       pow(gm,-0.25 );}
  else if (gm<2.0)  {return       pow(gm,-0.12 );}
  else              {return 0.0482025*pow((gm-2.0),2.0)+0.92019;}
}

double KikuchiAikenLRB_TCA::calcCHeqType1(double gm)
{
  if      (gm<2.0)  {return 1.0;}
  else              {return 0.036375*pow((gm-2.0),2.0)+1.0;}
}

double KikuchiAikenLRB_TCA::calcAType1(double gm, double heq, double u)
{
  if      (gm<2.0)  {return compABisection(heq,u,0.0,50.0,1e-6,26.501472);}
  else              {return 26.501472;}
}
// idac check print for OILES500R
//[Check print] hs_klrb
//  heq and u at 200% shear strain =   0.19093224874772838       0.31167639512587986
//  pa        at 200% shear strain =    26.501472000000042

double KikuchiAikenLRB_TCA::calcBType1(double gm, double a, double c, double heq, double u)
{
  if      (gm<2.0)  {return 0.0;}
  else              {return c*c*(M_PI*heq/u-(2+2/a*(exp(-2*a)-1)));}
}


//--------------------------
// type2 : OILES LRB250S
//--------------------------

double KikuchiAikenLRB_TCA::calcNType2(double gm)
{return 8.0;}

double KikuchiAikenLRB_TCA::calcPType2(double gm)
{
  if      (gm<2.0)  {return 0.0;}
  else              {return -0.35476 + 0.17738*gm;}
}

double KikuchiAikenLRB_TCA::calcCType2(double gm)
{return 6.0;}

double KikuchiAikenLRB_TCA::calcCQdType2(double gm)
{
  if      (gm<0.1)  {return 2.036*pow(gm,0.410);}
  else if (gm<0.5)  {return 1.106*pow(gm,0.145);}
  else              {return 1.0;}
}

double KikuchiAikenLRB_TCA::calcCKdType2(double gm)
{
  if      (gm<0.25) {return 0.779*pow(gm,-0.43 );}
  else if (gm<1.0)  {return       pow(gm,-0.25 );}
  else if (gm<2.0)  {return       pow(gm,-0.12 );}
  else              {return 0.0415275*pow((gm-2.0),2.0)+0.92019;}
}

double KikuchiAikenLRB_TCA::calcCHeqType2(double gm)
{
  if      (gm<2.0)  {return 1.0;}
  else              {return 0.0616*pow((gm-2.0),2.0)+1.0;}
}

double KikuchiAikenLRB_TCA::calcAType2(double gm, double heq, double u)
{
  if      (gm<2.0)  {return compABisection(heq,u,0.0,50.0,1e-6,26.501472);}
  else              {return 26.501472;}
}

// [Check print] hs_klrb for OILES250S ( this is the same valye as the one for OILES500R )
//  heq and u at 200% shear strain =   0.19093224874772838       0.31167639512587986
//  pa        at 200% shear strain =    26.501472000000042

double KikuchiAikenLRB_TCA::calcBType2(double gm, double a, double c, double heq, double u)
{
  if      (gm<2.0)  {return 0.0;}
  else              {return c*c*(M_PI*heq/u-(2+2/a*(exp(-2*a)-1)));}
}


//--------------------------
// type3 : Standard 400 (idac STANDARD-500-2�����j
//--------------------------

double KikuchiAikenLRB_TCA::calcNType3(double gm)
{return 8.0;}

double KikuchiAikenLRB_TCA::calcPType3(double gm)
{
  if      (gm<2.0)  {return 0.0;}
  else              {return -0.40074 + 0.20037*gm;}
}

double KikuchiAikenLRB_TCA::calcCType3(double gm)
{return 6.0;}

double KikuchiAikenLRB_TCA::calcCQdType3(double gm)
{
  if      (gm<0.1)  {return 2.036*pow(gm,0.410);}
  else if (gm<0.5)  {return 1.106*pow(gm,0.145);}
  else              {return 1.0;}
}

double KikuchiAikenLRB_TCA::calcCKdType3(double gm)
{
  if      (gm<0.25) {return 0.779*pow(gm,-0.43 );}
  else if (gm<1.0)  {return       pow(gm,-0.25 );}
  else if (gm<2.0)  {return       pow(gm,-0.12 );}
  else              {return 0.103415*pow((gm-2.0),2.0)+0.92019;}
}

double KikuchiAikenLRB_TCA::calcCHeqType3(double gm)
{
  if      (gm<2.0)  {return 1.0;}
  else              {return 0.036375*pow((gm-2.0),2.0)+1.0;}
}

double KikuchiAikenLRB_TCA::calcAType3(double gm, double heq, double u)
{
  if      (gm<2.0)  {return compABisection(heq,u,0.0,50.0,1e-6,26.501472);}
  else              {return 26.501472;}
}

// [Check print] hs_klrb for OILES250S ( this is the same valye as the one for OILES500R )
//  heq and u at 200% shear strain =   0.19093224874772838       0.31167639512587986
//  pa        at 200% shear strain =    26.501472000000042

double KikuchiAikenLRB_TCA::calcBType3(double gm, double a, double c, double heq, double u)
{
  if      (gm<2.0)  {return 0.0;}
  else              {return c*c*(M_PI*heq/u-(2+2/a*(exp(-2*a)-1)));}
}

