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

// Written:  Ippei Tsunezawa
// Created:  February 2014
// Modified: Feb 09, 2015
//
//
// Description: This file contains the class definition for EBSBearing.
//
// YamamotoBiaxialHDR.cppをひな形にして作成
//

#ifndef EBSBearing_h
#define EBSBearing_h

#include <Element.h>
#include <Matrix.h>

//#define ELE_TAG_EBSBearing 71816

class Channel;
class Response;

class EBSBearing : public Element
{
 public:
  // constructor
  EBSBearing(int Tag, int Nd1, int Nd2, int Tp, double DDo, double DDi, double Hr, double Mu, double Ec, double Gr, double Fz0, double As,
	  const Vector OriYp, const Vector OriX = 0,
	  double Mass = 0.0);

  EBSBearing();
  
  // destructor
  ~EBSBearing();
  
  // method to get class type
  const char *getClassType() const {return "EBSBearing";};
  
  // public methods to obtain information about dof & connectivity    
  int getNumExternalNodes() const;
  const ID &getExternalNodes();
  Node **getNodePtrs();
  int getNumDOF(); 
  void setDomain(Domain *theDomain);
  

  // public methods to set the state of the element    
  int commitState();
  int revertToLastCommit();        
  int revertToStart();      

  int update();
  
  // public methods to obtain stiffness, mass, damping and residual information    
  const Matrix &getTangentStiff();
  const Matrix &getInitialStiff();
  const Matrix &getMass();
  
  void zeroLoad();
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);
  
  const Vector &getResistingForce();
  const Vector &getResistingForceIncInertia();
  
  // public methods for element output
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  int displaySelf(Renderer &theViewer, int displayMode, float fact);    
  void Print(OPS_Stream &s, int flag = 0);    
  
  // public methods for element recorder
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInfo);
  
 protected:
  
 private:
  //-------------------------------------------------------------------------------
  int setTrialStrain(const Vector &strain);
  const double &getStrain(int direction);
  const double &getStress(int direction);
  const double &getTangent(int direction);
  const double &getInitialTangent(int direction);
  //-------------------------------------------------------------------------------

  // private methods
  void setUp();
  
  // private attributes - a copy for each object of the class
  ID connectedExternalNodes;        // contains the tags of the end nodes
  Node *theNodes[2];                // array of nodes
  //  UniaxialMaterial **theMaterials;  // materialをnSpring個

  // parameters
  //  int nSpring=1;
  Vector oriX;   // local x direction
  Vector oriYp;  // local yp direction
  double mass; // mass of element
  
  // 座標変換
  Matrix Tgl; // transformation matrix from global to local system
  Matrix Tlb; // transformation matrix from local to basic system
  
  // 変位・剛性
  Vector basicDisp;  //ub
  Vector localDisp;  //ul
  Vector basicForce; //qb
  Matrix basicStiff; //kb  
  Matrix basicStiffInit; //kbInit

  static Matrix theMatrix;
  static Vector theVector;
  static Vector theLoad;

  // 入力データ
  int tp; // モデルの種類
          // = 1; 弾性すべり支承
          // = 2; 弾性すべり支承(摩擦係数の面圧依存を考慮)

  double ddo; // 積層ゴム外径 [m]
  double ddi; // 積層ゴム内径 [m]
  double hr; // ゴム総厚 [m]
  double mu; //摩擦係数
  double ec; //圧縮弾性係数
  double gr; //せん断弾性係数
  double fz0;//初期軸力（圧縮方向は負）
  double as; //すべり材の断面積(m^2)

  // モデルパラメータ
  double ar;    // 積層ゴム断面積 [m^2]
  double DP[2];//変位増分
  double Kz;//鉛直ばね定数
  double Kxy;//せん断弾性率
  double Alpha;//滑り出し変位
  double Beta;//距離
  double P_AB;// Alpha/Beta
  double mui;//面圧に依存する摩擦係数


  // 初期剛性
  double initialStiff[3];

  // trial values
  double trialStiff[3];
  double trialDeform[3];
  double trialForce[3];
  double trialQ[3]; // 移動原点Qの中心

  // commit values
  double commitStiff[3];
  double commitDeform[3];
  double commitForce[3];
  double commitQ[3]; // 移動原点Qの中心
};

#endif
