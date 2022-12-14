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
// $Date: 2013-10-30 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Takeda.h,v $

// Written: Ippei Tsunezawa
// Created: Oct 30, 2013
// Modified: Feb 09, 2015
//
// Takeda model for Reinforced�]Concrete Materials
//
// Description: This file contains the class definition for Takeda.
//
//

#ifndef Takeda_h
#define Takeda_h

#include <UniaxialMaterial.h>

//#define MAT_TAG_Takeda 6109


class Takeda : public UniaxialMaterial
{
 public:
  Takeda(int tag, double fc, double fy, double ge, 
       double bc, double by, double alph);
  Takeda();
  ~Takeda();

  const char *getClassType(void) const {return "Takeda";};

  int setTrialStrain(double strain, double strainRate = 0.0);
  double getStrain(void);
  double getStress(void);
  double getTangent(void);
  double getInitialTangent(void);

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

  double trialDeform;  // trial deformation
  double trialForce;   // trial force
  double trialStiff;   // trial stiffness
  double commitDeform; // commit deformation
  double commitForce;  // commit force
  double commitStiff;  // commit stiffness

  //
  double fc;  // �N���b�N�_�׏d
  double fy;  // �~���_�׏d
  double ge;  // �e���o�l�̍���
  double bc;  // �����ቺ�W�� gc/ge
  double by;  // �����ቺ�W�� gy/ge
  double alph;// ���׎������ቺ�w��

  //
  double dm;  // (case3 unloading�ŗp����)
  double ga;  // (case4��case5�ŗp����)
  double gb;  // (case4��case5�ŗp����)
  double dum; // ���ד_Um�̕ψ�
  double fum; // ���ד_Um�̉׏d
  double gc;  // �N���b�N���̍���
  double gy;  // �~�����̍���
  double dc;  // �N���b�N�_�ψ�
  double dy;  // �~���_�ψ�
  double grv; // (fc+fy)/(dc+dy)

  // trial values
  int trialStg;     //���ɐi�ރX�e�[�W���L��
  int trialCrk;     //�N���b�N�̗L��
  double trialDmax; //�ψʂ̍ő�_
  double trialFmax; //�׏d�̍ő�_
  double trialDmin; //�ψʂ̍ŏ��_
  double trialFmin; //�׏d�̍ŏ��_
  double trialDu0;  //���ד_0(�ψ�)
  double trialFu0;  //���ד_0(�׏d)
  double trialDu1;  //���ד_1(�ψ�)
  double trialFu1;  //���ד_1(�׏d)
  double trialDu2;  //���ד_2(�ψ�)
  double trialFu2;  //���ד_2(�׏d)
  double trialDu3;  //���ד_3(�ψ�)
  double trialFu3;  //���ד_3(�׏d)
  double trialx0;   //�׏d�[���̓_0
  double trialx1;   //�׏d�[���̓_1
  double trialx2;   //�׏d�[���̓_2
  double trialx3;   //�׏d�[���̓_3
  double trialgss;  //unloading�̍���
  double trialDcPos;//(case2�ŗp����)
  double trialDcNeg;//(case2�ŗp����)
  double trialgc0;  //�N���b�N�_����~���_�̍���

  // commit values
  int commitStg;     //���ɐi�ރX�e�[�W���L��
  int commitCrk;     //�N���b�N�̗L��
  double commitDmax; //�ψʂ̍ő�_
  double commitFmax; //�׏d�̍ő�_
  double commitDmin; //�ψʂ̍ŏ��_
  double commitFmin; //�׏d�̍ŏ��_
  double commitDu0;  //���ד_0(�ψ�)
  double commitFu0;  //���ד_0(�׏d)
  double commitDu1;  //���ד_1(�ψ�)
  double commitFu1;  //���ד_1(�׏d)
  double commitDu2;  //���ד_2(�ψ�)
  double commitFu2;  //���ד_2(�׏d)
  double commitDu3;  //���ד_3(�ψ�)
  double commitFu3;  //���ד_3(�׏d)
  double commitx0;   //�׏d�[���̓_0
  double commitx1;   //�׏d�[���̓_1
  double commitx2;   //�׏d�[���̓_2
  double commitx3;   //�׏d�[���̓_3
  double commitgss;  //unloading�̍���
  double commitDcPos;//(case2�ŗp����)
  double commitDcNeg;//(case2�ŗp����)
  double commitgc0;  //�N���b�N�_����~���_�̍���
};

#endif
