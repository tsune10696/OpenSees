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
// $Date: 2013-05-31 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/MultipleShearSpring_forThermalConductivityAnalysis/MultipleShearSpring_forThermalConductivityAnalysis.h,v $


// Written: Ippei Tsunezawa
// Created: October 2014
//
// Multiple Shear Spring (MSS) model for Thermal Conductivity Analysis
//
// Description: This file contains the class definition for MultipleShearSpring_forThermalConductivityAnalysis.


#ifndef MultipleShearSpring_forThermalConductivityAnalysis_h
#define MultipleShearSpring_forThermalConductivityAnalysis_h

#include <Element.h>
#include <Matrix.h>


class Channel;
class KikuchiAikenLRB_TCA;
class Response;

class MultipleShearSpring_forThermalConductivityAnalysis : public Element
{
 public:
  // constructor
  MultipleShearSpring_forThermalConductivityAnalysis(int Tag, int Nd1, int Nd2,
		      int NSpring,
		      KikuchiAikenLRB_TCA *Material,
		      double LimDisp, double Dtstep,
		      int Tate, int Yoko, int Hyoko, int Ctate, int Ftate, int Jtate, int Itate,
		      int Norl, int Noss,
                      double Hc, double Hf, double Hj, double Hrl, double Hss, double Hi,
		      int Target_type,
		      int Targettate , int Targetyoko, int Loop, int Totalnumber,
                      double Capap, double Capaf, double Capac, double Caparl,
		      double Condp, double Condf, double Condc, double Condrl,
		      const Vector OriYp, const Vector OriX = 0,
		      double Mass = 0.0 );
  MultipleShearSpring_forThermalConductivityAnalysis();
  
  // destructor
  ~MultipleShearSpring_forThermalConductivityAnalysis();
  
  // method to get class type
  const char *getClassType() const {return "MultipleShearSpring_forThermalConductivityAnalysis";};
  
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

  //�ǉ���----------------
  double temp();
  int ThermalConductivityAnalysis();
  int SumEnergy();
  //---------------------
  
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
  // private methods
  void setUp();
  
  // private attributes - a copy for each object of the class
  ID connectedExternalNodes;        // contains the tags of the end nodes
  Node *theNodes[2];                // array of nodes
  KikuchiAikenLRB_TCA **theMaterials;  // materials
  
  // parameters
  int nSpring; //number of shear springs in MSS
  double *cosTht; //arrangement of each spring (cos)
  double *sinTht; //arrangement of each spring (sin)
  Vector oriX;   // local x direction
  Vector oriYp;  // local yp direction
  double mass; // mass of element
  Vector allTemperatures;                  //�S���x���



  //�M�`����͂ɗp����z��f�[�^���i�[����\����
  struct Data{
    double Trial_temp;                   //���̓_�̉��x
    double Commit_temp;                  //�O�̓_�̉��x
    double Trial_tempX;                  //���̓_�̉��x�i���̒l�j
    double Commit_tempX;                 //�O�̓_�̉��x�i���̒l�j
    double mesh_dr;                      //���̗v�f�̏��(mesh�̃�r[m])
    double mesh_dh;                      //���̗v�f�̏��(mesh�̃�h[m])
    int    mesh_mat;                     /*���̗v�f�̏��
					   (�������ʔԍ�(0= ��, 1= �V�R�S��, 2= �t�����W or �}���|�� or �A���|�� or �����|��, 3=�R���N���[�g))*/
    double heat_capacity;                //�M�e��
    double ruiseki_r;                    //�ݐ�r
    double ruiseki_h;                    //�ݐ�h
    double original_thermal_conductivity;//�����ŗL�̓����M�`����
    double volume;                       //�̐�
    double heat_capacity_by_volume;
    double boundary_right;               //�M�̍����P�ɂ܂Ƃ߂�
    double boundary_bottom;              //�M�̍����P�ɂ܂Ƃ߂�
  };

  Data **TCA_data;

  double center_distance_right;        //�Q�_�Ԃ̒��S�_�ԋ���
  double touch_area_right;             //�ڐG�ʐ�
  double thermal_conductivity_right;

  double center_distance_bottom;       //�Q�_�Ԃ̒��S�_�ԋ���
  double touch_area_bottom;            //�ڐG�ʐ�
  double thermal_conductivity_bottom;


  double *trialF;
  double *commitF;
  double *dDeform;
  double *trialD;
  double *commitD;


  int    tate;                                //�c�����i���������j�̕�����
  int    yoko;                                //�������i���a�����j�̕�����
  int    hyoko;                               //�������i���a�����j�̔��M�����̕�����
  int    ctate;                               //�c�����i���������j�̃R���N���[�g������
  int    ftate;                               //�c�����i���������j�̃t�����W������
  int    jtate;                               //�c�����i���������j�̘A���|������
  int    itate;                               //�c�����i���������j�̑}���|������

  int    norl;                                //�V�R�S���S�w��(�Ɛk���u�̎��ۂ̑���)
  int    noss;                                //�����|�S�w��(�Ɛk���u�̎��ۂ̑���)

  double hc;                                  //�R���N���[�g����[m]
  double hf;                                  //�t�����W����[m]
  double hj;                                  //�A���|����[m]
  double hrl;                                 //�V�R�S���P�w��[m]
  double hss;                                 //�����|�P�w��[m]
  double hi;                                  //�}���|����[m]

  double dtstep;                              //1�X�e�b�v�̎��ԊԊu

  int    targettate;                          //target_temp��x�Ԓn
  int    targetyoko;                          //target_temp��y�Ԓn

  int    loop;                                //1�X�e�b�v�ŉ񂵂��M�`����͂̉�

  int    totalnumber;                         //��͂ɗp���Ă���LRB�̈���(��͂̒P����)

  double capap;                               //���̔M�e��/�̐�[J/(m^3 * K)]
  double capaf;                               //�t�����W�̔M�e��/�̐�[J/(m^3 * K)]
  double capac;                               //�R���N���[�g�̔M�e��/�̐�[J/(m^3 * K)]
  double caparl;                              //�V�R�S���̔M�e��/�̐�[J/(m^3 * K)]
  double condp;                               //���̔M�`����[W/(m * K)]
  double condf;                               //�t�����W�̔M�`����[W/(m * K)]
  double condc;                               //�R���N���[�g�̔M�`����[W/(m * K)]
  double condrl;                              //�V�R�S���̔M�`����[W/(m * K)]

  int target_type;                            //�̗p���x�����i0=���b�V�����w�肵�Ă��̃��b�V���̉��x���펞�̗p, 1=���v���O�̍ő剷�x���̗p�j

  int    rltate;                              //�V�R�S���̏c�̗v�f���i�{���f���ł̗v�f���j
  int    sstate;                              //�����|�̏c�̗v�f���i�{���f���ł̗v�f���j


  // calculation of Feq and Seq
  double limDisp; //minimum deformation to calculate Feq and Seq (if limDisp is 0, never calculate)
  KikuchiAikenLRB_TCA *dmyMssMaterial; //imaginary material to calculate Feq and Seq
  double mssFeq; //equivalent coefficient for force
  double mssSeq; //equivalent coefficient for stiffness


  // transformation
  Matrix Tgl; // transformation matrix from global to local system
  Matrix Tlb; // transformation matrix from local to basic system
  
  // displacement, force, stiffness
  Vector basicDisp;  //ub
  Vector localDisp;  //ul
  Vector basicForce; //qb

  Vector basicQ2;

  Matrix basicStiff; //kb  
  Matrix basicStiffInit; //kbInit

  static Matrix theMatrix;
  static Vector theVector;
  static Vector theLoad;
  static Vector theHyst;

  double trialEsum;
  double commitEsum;
  double delta_Esum;

  double dt;                                  //�M�`����͂̎��ԊԊu
  double heated_part_ruiseki_r;               //�����b�V���̔��a���̑��a
  double heated_part_ruiseki_h;               //�����b�V���̍������̑��a
  double heated_part_r;                       //���M�����i���j�̔��a
  double heated_part_h;                       //���M�����i���j�̑�����
  double heated_part_volume;                  //���M�����i���j�̑̐�
  double zentai_r;                            //�ϑw�S���S�̔��a
  double lead_r;                              //���v���O���a
  double target_temp;
  double dfE;

};

#endif
