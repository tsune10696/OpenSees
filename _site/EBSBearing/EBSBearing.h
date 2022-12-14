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
// YamamotoBiaxialHDR.cpp���ЂȌ`�ɂ��č쐬
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
  //  UniaxialMaterial **theMaterials;  // material��nSpring��

  // parameters
  //  int nSpring=1;
  Vector oriX;   // local x direction
  Vector oriYp;  // local yp direction
  double mass; // mass of element
  
  // ���W�ϊ�
  Matrix Tgl; // transformation matrix from global to local system
  Matrix Tlb; // transformation matrix from local to basic system
  
  // �ψʁE����
  Vector basicDisp;  //ub
  Vector localDisp;  //ul
  Vector basicForce; //qb
  Matrix basicStiff; //kb  
  Matrix basicStiffInit; //kbInit

  static Matrix theMatrix;
  static Vector theVector;
  static Vector theLoad;

  // ���̓f�[�^
  int tp; // ���f���̎��
          // = 1; �e�����ׂ�x��
          // = 2; �e�����ׂ�x��(���C�W���̖ʈ��ˑ����l��)

  double ddo; // �ϑw�S���O�a [m]
  double ddi; // �ϑw�S�����a [m]
  double hr; // �S������ [m]
  double mu; //���C�W��
  double ec; //���k�e���W��
  double gr; //����f�e���W��
  double fz0;//�������́i���k�����͕��j
  double as; //���ׂ�ނ̒f�ʐ�(m^2)

  // ���f���p�����[�^
  double ar;    // �ϑw�S���f�ʐ� [m^2]
  double DP[2];//�ψʑ���
  double Kz;//�����΂˒萔
  double Kxy;//����f�e����
  double Alpha;//����o���ψ�
  double Beta;//����
  double P_AB;// Alpha/Beta
  double mui;//�ʈ��Ɉˑ����門�C�W��


  // ��������
  double initialStiff[3];

  // trial values
  double trialStiff[3];
  double trialDeform[3];
  double trialForce[3];
  double trialQ[3]; // �ړ����_Q�̒��S

  // commit values
  double commitStiff[3];
  double commitDeform[3];
  double commitForce[3];
  double commitQ[3]; // �ړ����_Q�̒��S
};

#endif
