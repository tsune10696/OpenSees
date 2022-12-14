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
// $Date: 2012-08-09 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/KikuchiAikenLRB_TCA.h,v $

// Written: Ippei Tsunezawa
// Created: February 1, 2015
//
// Kikuchi&Aiken model for lead rubber bearing (for TCA)
//
// Description: This file contains the class definition for KikuchiAikenLRB_TCA.
//
// �����𑥂����̎菇�Ōv�Z����
// (����)                        (�o��)
//  �ψ� -> �Ђ��� ->  ����  ->   �׏d
//          �Ђ��� -> �e���� -> �΂˒萔
//


#ifndef KikuchiAikenLRB_TCA_h
#define KikuchiAikenLRB_TCA_h

#include <UniaxialMaterial.h>


class KikuchiAikenLRB_TCA : public UniaxialMaterial
{
 public:
  KikuchiAikenLRB_TCA(int tag, int type, double ar, double hr, double gr, double ap, double tp,
		  double alph, double beta, double temp, double rk, double rq, double rs, double rf);

  KikuchiAikenLRB_TCA();
  ~KikuchiAikenLRB_TCA();

  const char *getClassType(void) const {return "KikuchiAikenLRB_TCA";};

  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void);

//--------------------
  double getQ1(void);
  double getQ2(void);
  double getHr(void);
  double getAr(void);
  double getAp(void);
  double getQd(void);
  int setTemp(double target_temp);
  double getTemp(void);
//--------------------

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  UniaxialMaterial *getCopy(void);

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel,
	       FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag =0);

 protected:

 private:

  // ���̓f�[�^�ϐ�
  int    Type; // LRB�������̃^�C�v
               // =1: LRB500R (idac: OILES500F�����j
               // =2: LRB250S (idac: OILES250S�����j
               // =3: Standard400 (idac: STANDARD-500-2�����j
  double Ar;   // �ϑw�S���̒f�ʐ� [m^2]
  double Hr;   // �ϑw�S���̃S������ [m]
  double Gr;   // �S���̂���f�e���� [N/m^2]
  double Ap;   // ���v���O�̒f�ʐ� [m^2]
  double Tp;   // ���v���O�̍~�����͓x [N/m^2]
  double Alph; // ���v���O�݂̂����̂���f�e���� [N/m^2]
  double Beta; // ���������̍~���㍄���ɑ΂���{��
  double Temp; // ���x [��]
  double Rk;   // Kd�̕ϓ���
  double Rq;   // Qd�̕ϓ���
  double Rs;   // �����̒ጸ�� for MSS model
  double Rf;   // �׏d�̒ጸ�� for MSS model


  // �����͓����l
  double qd100; // �ؕЉ׏d�̊�l
  double kd100; // �~���㍄���̊�l
  double ku100; // ���������̊�l
  double qd;    // �ؕЉ׏d
  double kd;    // �~���㍄��
  double ku;    // ��������

  // ���E�Ђ��݁E��������
  double trgStrain; // �e�����E�Ђ���
  double lmtStrain; // �K�p���E�Ђ���
  double initialStiff;   // �ό`���Ȃ��Ƃ��̂΂˒萔

  // temporary values: �e�n���f���p�����[�^
  double tmpStrain; //�p�����[�^�Z�o�p�Ђ���
  double tmpDeform; //�p�����[�^�Z�o�p�ό`
  double keq; //keq
  double heq; //heq
  double u;   //u
  double n;   //n
  double p;   //p
  double a;   //a
  double b;   //b
  double c;   //c
  double xm;  //�ő�ψ�
  double fm;  //�ő�׏d
  double x;   //����ψ�
  double alpha; //�C���W��
  //double q1Stf; //Q1�̌X��
  //double q2Stf; //Q2�̌X��

  //-------------
  double Alph_T;
  double Tp_TCA;
  double Alph_T_ini;
  double Tp_ini;
  double ce;
  double trialQ1Stf; //Q1�̌X��
  double trialQ2Stf; //Q2�̌X��
  double commitQ1Stf; //Q1�̌X��
  double commitQ2Stf; //Q2�̌X��
  //-------------

  // trial values
  double trialDeform;     //���͗p, �ψ�
  double trialForce;           //�o�͗p, �׏d
  double trialStiff;       //�o�͗p, �΂˒萔
  double trialStrain;          //�����v�Z�p, �Ђ���
  bool   trialIfElastic;       //���`���E�t���O
  double trialQ1;              //����`�e������(����)
  double trialQ2;              //������������(����)
  double trialMaxStrain;       //�o���ő�Ђ���(��Βl)
  double trialDDeform;         //�ψʑ���
  int    trialDDeformLastSign; //�O��̕ψʑ����̕���,�ړ��̂Ȃ��ꍇ�͍X�V���Ȃ�
  int    trialIdxRev;          //�ŐV�̗����Ȑ��C���f�b�N�X


  // commit values
  double commitDeform;     //���͗p, �ψ�
  double commitForce;           //�o�͗p, �׏d
  double commitStiff;       //�o�͗p, �΂˒萔
  double commitStrain;          //�����v�Z�p, �Ђ���
  bool   commitIfElastic;       // ���`���E�t���O
  double commitQ1;              //����`�e������(����)
  double commitQ2;              //������������(����)
  double commitMaxStrain;       //�o���ő�Ђ���(��Βl)
  double commitDDeform;         //�Ђ��ݑ���
  int    commitDDeformLastSign; //�O��̂Ђ��ݑ����̕���,�ړ��̂Ȃ��ꍇ�͍X�V���Ȃ�
  int    commitIdxRev;          //�ŐV�̗����Ȑ��C���f�b�N�X


  // �����Ȑ��L���p�ϐ�
  // XX[0] �X�P���g���J�[�u
  // XX[1] �X�P���g���J�[�u����̏���
  // XX[2,3,...] ���]��
  int numIdx; //�L�������:����500��,�ǉ�500����
  double *revXBgn;  //���]�J�n�_��x
  double *revQ2Bgn; //���]�J�n�_��q2
  double *revXEnd;  //���]�w���_��x
  double *revQ2End; //���]�w���_��q2
  double *revB;     //���ׁE���]����b
  double *revAlpha; //���]���̏C���W��


  //�e�n���f���Z�莮(�S����ɂ�炸����)
  double compQ1(double u, double n, double p, double fm, double x);
  double compQ2Unload(double u, double a, double b, double c, double fm, double x);
  double compQ2Masing(double u, double a, double b, double c, double fm, double x1, double x2, double q2i, double alpha);
  double compAlpha(double a, double b1, double b2, double c, double x1, double x2, double alpha0);

  double compQ1Derivertive(double u, double n, double p, double keq, double x);
  double compQ2UnloadDerivertive(double u, double a, double b, double c, double keq, double x);
  double compQ2MasingDerivertive(double u, double a, double b, double c, double keq, double x1, double x2, double alpha);

  //�p�����[�^a�Z�o�p�̓񕪖@
  static double compABisection(double heq, double u, double min, double max, double lim, double tol);

  //�����A�����萔�̌v�Z
  static double compKeq(double xm, double qd, double kd);
  static double compHeq(double xm, double qd, double kd, double ku);

  //�e�p�����[�^�̎Z�莮(�S���킲�Ƃɒ�`�������̂��|�C���^�őI��)
  double (*calcN)(double gm);
  double (*calcP)(double gm);
  double (*calcA)(double gm, double heq, double u);
  double (*calcB)(double gm, double a, double c,double heq, double u);
  double (*calcC)(double gm);
  double (*calcCQd)(double gm);
  double (*calcCKd)(double gm);
  double (*calcCHeq)(double gm);


  //�S����ɉ������e��֐��̒�`
  //OILES LRB500R
  static double calcNType1(double gm);
  static double calcPType1(double gm);
  static double calcAType1(double gm, double heq, double u);
  static double calcBType1(double gm, double a, double c,double heq, double u);
  static double calcCType1(double gm);
  static double calcCQdType1(double gm);
  static double calcCKdType1(double gm);
  static double calcCHeqType1(double gm);
  //OILES LRB250S
  static double calcNType2(double gm);
  static double calcPType2(double gm);
  static double calcAType2(double gm, double heq, double u);
  static double calcBType2(double gm, double a, double c,double heq, double u);
  static double calcCType2(double gm);
  static double calcCQdType2(double gm);
  static double calcCKdType2(double gm);
  static double calcCHeqType2(double gm);
  //STANDARD 400
  static double calcNType3(double gm);
  static double calcPType3(double gm);
  static double calcAType3(double gm, double heq, double u);
  static double calcBType3(double gm, double a, double c,double heq, double u);
  static double calcCType3(double gm);
  static double calcCQdType3(double gm);
  static double calcCKdType3(double gm);
  static double calcCHeqType3(double gm);

};

#endif

