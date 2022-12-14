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
// YamamotoBiaxialHDR.cpp‚ğ‚Ğ‚ÈŒ`‚É‚µ‚Äì¬
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
  //  UniaxialMaterial **theMaterials;  // material‚ğnSpringŒÂ

  // parameters
  //  int nSpring=1;
  Vector oriX;   // local x direction
  Vector oriYp;  // local yp direction
  double mass; // mass of element
  
  // À•W•ÏŠ·
  Matrix Tgl; // transformation matrix from global to local system
  Matrix Tlb; // transformation matrix from local to basic system
  
  // •ÏˆÊE„«
  Vector basicDisp;  //ub
  Vector localDisp;  //ul
  Vector basicForce; //qb
  Matrix basicStiff; //kb  
  Matrix basicStiffInit; //kbInit

  static Matrix theMatrix;
  static Vector theVector;
  static Vector theLoad;

  // “ü—Íƒf[ƒ^
  int tp; // ƒ‚ƒfƒ‹‚Ìí—Ş
          // = 1; ’e«‚·‚×‚èx³
          // = 2; ’e«‚·‚×‚èx³(–€CŒW”‚Ì–Êˆ³ˆË‘¶‚ğl—¶)

  double ddo; // Ï‘wƒSƒ€ŠOŒa [m]
  double ddi; // Ï‘wƒSƒ€“àŒa [m]
  double hr; // ƒSƒ€‘Œú [m]
  double mu; //–€CŒW”
  double ec; //ˆ³k’e«ŒW”
  double gr; //‚¹‚ñ’f’e«ŒW”
  double fz0;//‰Šú²—Íiˆ³k•ûŒü‚Í•‰j
  double as; //‚·‚×‚èŞ‚Ì’f–ÊÏ(m^2)

  // ƒ‚ƒfƒ‹ƒpƒ‰ƒ[ƒ^
  double ar;    // Ï‘wƒSƒ€’f–ÊÏ [m^2]
  double DP[2];//•ÏˆÊ‘•ª
  double Kz;//‰”’¼‚Î‚Ë’è”
  double Kxy;//‚¹‚ñ’f’e«—¦
  double Alpha;//ŠŠ‚èo‚µ•ÏˆÊ
  double Beta;//‹——£
  double P_AB;// Alpha/Beta
  double mui;//–Êˆ³‚ÉˆË‘¶‚·‚é–€CŒW”


  // ‰Šú„«
  double initialStiff[3];

  // trial values
  double trialStiff[3];
  double trialDeform[3];
  double trialForce[3];
  double trialQ[3]; // ˆÚ“®Œ´“_Q‚Ì’†S

  // commit values
  double commitStiff[3];
  double commitDeform[3];
  double commitForce[3];
  double commitQ[3]; // ˆÚ“®Œ´“_Q‚Ì’†S
};

#endif
