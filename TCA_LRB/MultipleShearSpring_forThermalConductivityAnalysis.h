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

  //追加分----------------
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
  Vector allTemperatures;                  //全温度情報



  //熱伝導解析に用いる配列データを格納する構造体
  struct Data{
    double Trial_temp;                   //今の点の温度
    double Commit_temp;                  //前の点の温度
    double Trial_tempX;                  //今の点の温度（仮の値）
    double Commit_tempX;                 //前の点の温度（仮の値）
    double mesh_dr;                      //その要素の情報(meshのΔr[m])
    double mesh_dh;                      //その要素の情報(meshのΔh[m])
    int    mesh_mat;                     /*その要素の情報
					   (物質識別番号(0= 鉛, 1= 天然ゴム, 2= フランジ or 挿入鋼板 or 連結鋼板 or 内部鋼板, 3=コンクリート))*/
    double heat_capacity;                //熱容量
    double ruiseki_r;                    //累積r
    double ruiseki_h;                    //累積h
    double original_thermal_conductivity;//物質固有の等価熱伝導率
    double volume;                       //体積
    double heat_capacity_by_volume;
    double boundary_right;               //熱の項を１つにまとめた
    double boundary_bottom;              //熱の項を１つにまとめた
  };

  Data **TCA_data;

  double center_distance_right;        //２点間の中心点間距離
  double touch_area_right;             //接触面積
  double thermal_conductivity_right;

  double center_distance_bottom;       //２点間の中心点間距離
  double touch_area_bottom;            //接触面積
  double thermal_conductivity_bottom;


  double *trialF;
  double *commitF;
  double *dDeform;
  double *trialD;
  double *commitD;


  int    tate;                                //縦方向（高さ方向）の分割数
  int    yoko;                                //横方向（半径方向）の分割数
  int    hyoko;                               //横方向（半径方向）の発熱部分の分割数
  int    ctate;                               //縦方向（高さ方向）のコンクリート分割数
  int    ftate;                               //縦方向（高さ方向）のフランジ分割数
  int    jtate;                               //縦方向（高さ方向）の連結鋼板分割数
  int    itate;                               //縦方向（高さ方向）の挿入鋼板分割数

  int    norl;                                //天然ゴム全層数(免震装置の実際の総数)
  int    noss;                                //内部鋼板全層数(免震装置の実際の総数)

  double hc;                                  //コンクリート高さ[m]
  double hf;                                  //フランジ高さ[m]
  double hj;                                  //連結鋼板高さ[m]
  double hrl;                                 //天然ゴム１層厚[m]
  double hss;                                 //内部鋼板１層厚[m]
  double hi;                                  //挿入鋼板高さ[m]

  double dtstep;                              //1ステップの時間間隔

  int    targettate;                          //target_tempのx番地
  int    targetyoko;                          //target_tempのy番地

  int    loop;                                //1ステップで回した熱伝導解析の回数

  int    totalnumber;                         //解析に用いているLRBの員数(解析の単純化)

  double capap;                               //鉛の熱容量/体積[J/(m^3 * K)]
  double capaf;                               //フランジの熱容量/体積[J/(m^3 * K)]
  double capac;                               //コンクリートの熱容量/体積[J/(m^3 * K)]
  double caparl;                              //天然ゴムの熱容量/体積[J/(m^3 * K)]
  double condp;                               //鉛の熱伝導率[W/(m * K)]
  double condf;                               //フランジの熱伝導率[W/(m * K)]
  double condc;                               //コンクリートの熱伝導率[W/(m * K)]
  double condrl;                              //天然ゴムの熱伝導率[W/(m * K)]

  int target_type;                            //採用温度取り方（0=メッシュを指定してそのメッシュの温度を常時採用, 1=鉛プラグの最大温度を採用）

  int    rltate;                              //天然ゴムの縦の要素数（本モデルでの要素数）
  int    sstate;                              //内部鋼板の縦の要素数（本モデルでの要素数）


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

  double dt;                                  //熱伝導解析の時間間隔
  double heated_part_ruiseki_r;               //鉛メッシュの半径情報の総和
  double heated_part_ruiseki_h;               //鉛メッシュの高さ情報の総和
  double heated_part_r;                       //発熱部分（鉛）の半径
  double heated_part_h;                       //発熱部分（鉛）の総高さ
  double heated_part_volume;                  //発熱部分（鉛）の体積
  double zentai_r;                            //積層ゴム全体半径
  double lead_r;                              //鉛プラグ半径
  double target_temp;
  double dfE;

};

#endif
