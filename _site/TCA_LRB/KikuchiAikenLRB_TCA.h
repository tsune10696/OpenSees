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
// ※履歴則を次の手順で計算する
// (入力)                        (出力)
//  変位 -> ひずみ ->  応力  ->   荷重
//          ひずみ -> 弾性率 -> ばね定数
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

  // 入力データ変数
  int    Type; // LRB実験式のタイプ
               // =1: LRB500R (idac: OILES500F相当）
               // =2: LRB250S (idac: OILES250S相当）
               // =3: Standard400 (idac: STANDARD-500-2相当）
  double Ar;   // 積層ゴムの断面積 [m^2]
  double Hr;   // 積層ゴムのゴム総厚 [m]
  double Gr;   // ゴムのせん断弾性率 [N/m^2]
  double Ap;   // 鉛プラグの断面積 [m^2]
  double Tp;   // 鉛プラグの降伏応力度 [N/m^2]
  double Alph; // 鉛プラグのみかけのせん断弾性率 [N/m^2]
  double Beta; // 初期剛性の降伏後剛性に対する倍率
  double Temp; // 温度 [℃]
  double Rk;   // Kdの変動率
  double Rq;   // Qdの変動率
  double Rs;   // 剛性の低減率 for MSS model
  double Rf;   // 荷重の低減率 for MSS model


  // 復元力特性値
  double qd100; // 切片荷重の基準値
  double kd100; // 降伏後剛性の基準値
  double ku100; // 初期剛性の基準値
  double qd;    // 切片荷重
  double kd;    // 降伏後剛性
  double ku;    // 初期剛性

  // 限界ひずみ・初期剛性
  double trgStrain; // 弾性限界ひずみ
  double lmtStrain; // 適用限界ひずみ
  double initialStiff;   // 変形がないときのばね定数

  // temporary values: 菊地モデルパラメータ
  double tmpStrain; //パラメータ算出用ひずみ
  double tmpDeform; //パラメータ算出用変形
  double keq; //keq
  double heq; //heq
  double u;   //u
  double n;   //n
  double p;   //p
  double a;   //a
  double b;   //b
  double c;   //c
  double xm;  //最大変位
  double fm;  //最大荷重
  double x;   //基準化変位
  double alpha; //修正係数
  //double q1Stf; //Q1の傾き
  //double q2Stf; //Q2の傾き

  //-------------
  double Alph_T;
  double Tp_TCA;
  double Alph_T_ini;
  double Tp_ini;
  double ce;
  double trialQ1Stf; //Q1の傾き
  double trialQ2Stf; //Q2の傾き
  double commitQ1Stf; //Q1の傾き
  double commitQ2Stf; //Q2の傾き
  //-------------

  // trial values
  double trialDeform;     //入力用, 変位
  double trialForce;           //出力用, 荷重
  double trialStiff;       //出力用, ばね定数
  double trialStrain;          //内部計算用, ひずみ
  bool   trialIfElastic;       //線形限界フラグ
  double trialQ1;              //非線形弾性成分(応力)
  double trialQ2;              //履歴減衰成分(応力)
  double trialMaxStrain;       //経験最大ひずみ(絶対値)
  double trialDDeform;         //変位増分
  int    trialDDeformLastSign; //前回の変位増分の符号,移動のない場合は更新しない
  int    trialIdxRev;          //最新の履歴曲線インデックス


  // commit values
  double commitDeform;     //入力用, 変位
  double commitForce;           //出力用, 荷重
  double commitStiff;       //出力用, ばね定数
  double commitStrain;          //内部計算用, ひずみ
  bool   commitIfElastic;       // 線形限界フラグ
  double commitQ1;              //非線形弾性成分(応力)
  double commitQ2;              //履歴減衰成分(応力)
  double commitMaxStrain;       //経験最大ひずみ(絶対値)
  double commitDDeform;         //ひずみ増分
  int    commitDDeformLastSign; //前回のひずみ増分の符号,移動のない場合は更新しない
  int    commitIdxRev;          //最新の履歴曲線インデックス


  // 履歴曲線記憶用変数
  // XX[0] スケルトンカーブ
  // XX[1] スケルトンカーブからの除荷
  // XX[2,3,...] 反転中
  int numIdx; //記憶する個数:初期500個,追加500個ずつ
  double *revXBgn;  //反転開始点のx
  double *revQ2Bgn; //反転開始点のq2
  double *revXEnd;  //反転指向点のx
  double *revQ2End; //反転指向点のq2
  double *revB;     //除荷・反転時のb
  double *revAlpha; //反転時の修正係数


  //菊地モデル算定式(ゴム種によらず共通)
  double compQ1(double u, double n, double p, double fm, double x);
  double compQ2Unload(double u, double a, double b, double c, double fm, double x);
  double compQ2Masing(double u, double a, double b, double c, double fm, double x1, double x2, double q2i, double alpha);
  double compAlpha(double a, double b1, double b2, double c, double x1, double x2, double alpha0);

  double compQ1Derivertive(double u, double n, double p, double keq, double x);
  double compQ2UnloadDerivertive(double u, double a, double b, double c, double keq, double x);
  double compQ2MasingDerivertive(double u, double a, double b, double c, double keq, double x1, double x2, double alpha);

  //パラメータa算出用の二分法
  static double compABisection(double heq, double u, double min, double max, double lim, double tol);

  //剛性、減衰定数の計算
  static double compKeq(double xm, double qd, double kd);
  static double compHeq(double xm, double qd, double kd, double ku);

  //各パラメータの算定式(ゴム種ごとに定義したものをポインタで選択)
  double (*calcN)(double gm);
  double (*calcP)(double gm);
  double (*calcA)(double gm, double heq, double u);
  double (*calcB)(double gm, double a, double c,double heq, double u);
  double (*calcC)(double gm);
  double (*calcCQd)(double gm);
  double (*calcCKd)(double gm);
  double (*calcCHeq)(double gm);


  //ゴム種に応じた各種関数の定義
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

