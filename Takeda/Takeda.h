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
// Takeda model for Reinforced‐Concrete Materials
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
  double fc;  // クラック点荷重
  double fy;  // 降伏点荷重
  double ge;  // 弾性バネの剛性
  double bc;  // 剛性低下係数 gc/ge
  double by;  // 剛性低下係数 gy/ge
  double alph;// 除荷時剛性低下指数

  //
  double dm;  // (case3 unloadingで用いる)
  double ga;  // (case4とcase5で用いる)
  double gb;  // (case4とcase5で用いる)
  double dum; // 除荷点Umの変位
  double fum; // 除荷点Umの荷重
  double gc;  // クラック時の剛性
  double gy;  // 降伏時の剛性
  double dc;  // クラック点変位
  double dy;  // 降伏点変位
  double grv; // (fc+fy)/(dc+dy)

  // trial values
  int trialStg;     //次に進むステージを記憶
  int trialCrk;     //クラックの有無
  double trialDmax; //変位の最大点
  double trialFmax; //荷重の最大点
  double trialDmin; //変位の最小点
  double trialFmin; //荷重の最小点
  double trialDu0;  //除荷点0(変位)
  double trialFu0;  //除荷点0(荷重)
  double trialDu1;  //除荷点1(変位)
  double trialFu1;  //除荷点1(荷重)
  double trialDu2;  //除荷点2(変位)
  double trialFu2;  //除荷点2(荷重)
  double trialDu3;  //除荷点3(変位)
  double trialFu3;  //除荷点3(荷重)
  double trialx0;   //荷重ゼロの点0
  double trialx1;   //荷重ゼロの点1
  double trialx2;   //荷重ゼロの点2
  double trialx3;   //荷重ゼロの点3
  double trialgss;  //unloadingの剛性
  double trialDcPos;//(case2で用いる)
  double trialDcNeg;//(case2で用いる)
  double trialgc0;  //クラック点から降伏点の剛性

  // commit values
  int commitStg;     //次に進むステージを記憶
  int commitCrk;     //クラックの有無
  double commitDmax; //変位の最大点
  double commitFmax; //荷重の最大点
  double commitDmin; //変位の最小点
  double commitFmin; //荷重の最小点
  double commitDu0;  //除荷点0(変位)
  double commitFu0;  //除荷点0(荷重)
  double commitDu1;  //除荷点1(変位)
  double commitFu1;  //除荷点1(荷重)
  double commitDu2;  //除荷点2(変位)
  double commitFu2;  //除荷点2(荷重)
  double commitDu3;  //除荷点3(変位)
  double commitFu3;  //除荷点3(荷重)
  double commitx0;   //荷重ゼロの点0
  double commitx1;   //荷重ゼロの点1
  double commitx2;   //荷重ゼロの点2
  double commitx3;   //荷重ゼロの点3
  double commitgss;  //unloadingの剛性
  double commitDcPos;//(case2で用いる)
  double commitDcNeg;//(case2で用いる)
  double commitgc0;  //クラック点から降伏点の剛性
};

#endif
