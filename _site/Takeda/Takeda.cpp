/* ****************************************************************** ***
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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Takeda.cpp,v $

// Written: Ippei Tsunezawa
// Created: Oct 30, 2013
// Modified: Feb 09, 2015
//
// Takeda model for Reinforced‐Concrete Materials
//
// Description: This file contains the class definition for Takeda and  the function to parse the TCL input.
// uniaxialMaterial Takeda matTag? fc? fy? ge? bc? by? alph?

#include <TclModelBuilder.h>
#include <Takeda.h>

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <tcl.h>


int
TclCommand_Takeda(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  //最低限必要な引数
  int tag;
  double fc;
  double fy;
  double ge;
  double bc;
  double by;
  double alph;

      
  //作成するオブジェクト
  UniaxialMaterial *theMaterial = 0;


  //入力引数のエラーチェック
  bool ifNoError = true;
  
  if (argc < 9) { // uniaxialMaterial Takeda matTag? fc? fy? ge? bc? by? alph?

    opserr << "WARNING Takeda invalid number of arguments\n";
    ifNoError = false;
  }


    //argv[2~5]のチェック
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING Takeda invalid tag" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[3], &fc ) != TCL_OK || fc <= 0.0 ) {
      opserr << "WARNING Takeda invalid fc" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[4], &fy) != TCL_OK || fy <= 0.0 ) {
      opserr << "WARNING Takeda invalid fy" << endln;
      ifNoError = false;
    }

    if (fy < fc) {
      opserr << "WARNING Takeda invalid fy and fc" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[5], &ge) != TCL_OK) {
      opserr << "WARNING Takeda invalid ge" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[6], &bc) != TCL_OK || bc <= 0.0 ) {
      opserr << "WARNING Takeda invalid bc" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[7], &by) != TCL_OK || by <= 0.0 ) {
      opserr << "WARNING Takeda invalid by" << endln;
      ifNoError = false;
    }

    if (by > bc) {
      opserr << "WARNING Takeda invalid by and bc" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[8], &alph) != TCL_OK || alph < 0.0 ) {
      opserr << "WARNING Takeda invalid alph" << endln;
      ifNoError = false;
    }



  //入力引数のエラー処理
  if (!ifNoError) {
    //入力データ
    opserr << "Input command: ";
    for (int i=0; i<argc; i++){
      opserr << argv[i] << " ";
    }
    opserr << endln;
    
    //必要データ
    opserr << "Want: uniaxialMaterial Takeda matTag? fc? fy? ge? bc? by? alph?" << endln;
    return TCL_ERROR;
  }


  //デフォルト値を設定しようとして間違えた場合に対応
  if (alph == 0.0) alph = 0.4;


  
//   // check print-----------------------------------
//   opserr << "Test1\n";
//   opserr << " argc= " << argc << endln;
//   opserr << " argv[1]= " << argv[1] <<endln;
//   opserr << " argv[2]= " << argv[2] << endln;
//   opserr << " argv[3]= " << argv[3] << endln;
//   opserr << " argv[4]= " << argv[4] << endln;
//   opserr << " argv[5]= " << argv[5] << endln;
//   opserr << " argv[6]= " << argv[6] << endln;
//   opserr << " argv[7]= " << argv[7] << endln;
//   opserr << " argv[8]= " << argv[8] << endln;
//   // ----------------------------------------------

//   // check print-----------------------------------
//   opserr << "Test\n";
//   opserr << " argc= " << argc << endln;
//   opserr << " tag= " << tag <<endln;
//   opserr << " fc= " << fc <<endln;
//   opserr << " fy= " << fy << endln;
//   opserr << " ge= " << ge << endln;
//   opserr << " bc= " << bc << endln;
//   opserr << " by= " << by << endln;
//   opserr << " alph= " << alph << endln;
//   // ----------------------------------------------
 


  
  // Parsing was successful, allocate the material
  theMaterial = new Takeda(tag, fc, fy, ge, bc, by, alph);

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




Takeda::Takeda(int tag, double fc, double fy, double ge, 
       double bc, double by, double alph)
  :UniaxialMaterial(tag,MAT_TAG_Takeda),fc(fc),fy(fy),ge(ge),
   bc(bc),by(by),alph(alph)
{

/////////////////////////////////////////////////////////////////
//  opserr << "  ---Takeda::Takeda  " << endln; 
/////////////////////////////////////////////////////////////////

//revertToStartの値を用いて、制御データを初期化する

revertToStart();

}


Takeda::Takeda()
  :UniaxialMaterial(0,MAT_TAG_Takeda)
{
  trialDeform  = 0.0;
  trialForce   = 0.0;
  trialStiff   = ge;
  commitDeform = 0.0;
  commitForce  = 0.0;
  commitStiff  = ge;
}

Takeda::~Takeda()
{

}

int
Takeda::setTrialStrain(double strain, double strainRate)
{

  // 履歴則  変位から荷重を計算する

  trialDeform = strain;
  double dDeform = trialDeform - commitDeform;


  // 仮荷重
  trialForce = commitForce + commitStiff * dDeform;


   // check print--------------------------------------------
   //   opserr << "   \n";
   //   opserr << "---Takeda::setTrialStrain start\n";
   //   opserr << "  trialDeform = " << trialDeform << endln;
   //   opserr << "  commitDeform= " << commitDeform << endln;
   //   opserr << "  dDeform     = " << dDeform << endln;

   //  opserr << "  trialForce     = " << trialForce  << endln;
   //  opserr << "  trialDeform    = " << trialDeform << endln;
   //  opserr << "  trialStiff     = " << trialStiff  << endln;
   //  opserr << "  trialStg       = " << trialStg    << endln;
   //  opserr << "  trialCrk       = " << trialCrk    << endln;
   //  opserr << "  commitStg      = " << commitStg   << endln;

   // -------------------------------------------------------


  switch(commitStg){

  //弾性状態
  case 1:

    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {

      //荷重がクラック点に達していない
      if(fabs(trialForce) <= fc) {
        trialStiff = ge;
        trialForce = ge * trialDeform;
        trialStg = 1;

      }

      //荷重がクラック点を超える
      else {


        //荷重が降伏点に達していない
        if(fabs(trialForce) <= fy){
          if(trialForce > 0.0) {
            trialStiff = commitgc0;
            trialForce = fc + commitgc0 * (trialDeform - commitDcPos); 
            trialStg = 2;
            trialCrk = commitCrk + 1;
          }
          else {
            trialStiff = commitgc0;
            trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
            trialStg = 2;
            trialCrk = commitCrk + 1;
          }
        }
        //荷重が降伏点を超える
        else {
          if(trialForce > 0.0) {
            trialStiff = gy;
            trialForce = fy + gy * (trialDeform - dy);
            trialStg = 3;
            trialCrk = commitCrk + 1;
          }
          else {
            trialStiff = gy;
            trialForce = -fy + gy * (trialDeform + dy);
            trialStg = 3;
            trialCrk = commitCrk + 1;
          }
        }
      }
    }

    else {
    trialStiff = ge;
    trialForce = ge * trialDeform;
    trialStg = 1;
    }

    break;


  //クラックから降伏の初期カーブ上
  case 2:


    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {
        
      //荷重が降伏点に達していない
      if(fabs(trialForce) <= fy){
        if(trialForce > 0.0) {
          trialStiff = commitgc0;
          trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
          trialStg = 2;
        }
        else {
          trialStiff = commitgc0;
          trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
          trialStg = 2; 
        }
      }

      //荷重が降伏点を超える
      else {
        if(trialForce > 0.0) {
          trialStiff = gy;
          trialForce = fy + gy * (trialDeform - dy);
          trialStg = 3;
        }
        else {
          trialStiff = gy;
            trialForce = -fy + gy * (trialDeform + dy);
            trialStg = 3;
        }
      }
    }
    else{
      if(trialForce > 0.0){
        if(trialDmax <= commitDeform)
          trialDmax = commitDeform;
        if(trialFmax <= commitForce)
          trialFmax = commitForce;

        trialStiff = (commitForce + fc) / (commitDeform + dc);
        trialForce = commitForce + trialStiff * (trialDeform - commitDeform);
        trialgss = trialStiff;
        trialStg = 5;
      }
      else {
        if(trialDmin >= commitDeform)
          trialDmin = commitDeform;
        if(trialFmin >= commitForce)
          trialFmin = commitForce;
        trialStiff = (-commitForce + fc) / (-commitDeform + dc);
        trialForce = commitForce + trialStiff * (trialDeform - commitDeform);
        trialgss = trialStiff;
        trialStg = 5;
      }
    }
    break;




  //降伏後の初期カーブ上の載荷
  case 3:

    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {
      if(trialForce > 0.0) {
        trialStiff = gy;
        trialForce = fy + gy * (trialDeform - dy);
        trialStg = 3;
      }
      else {
        trialStiff = gy;
        trialForce = -fy + gy * (trialDeform + dy);
        trialStg = 3;
      }
    }


    //unloading
    else {
      if(trialForce > 0.0) {
        if(trialDmax <= commitDeform)
          trialDmax = commitDeform;
        if(trialFmax <= commitForce)
          trialFmax = commitForce;
        dm = fabs(trialDmax);
        if(fabs(trialDmin) > dm)
          dm = fabs(trialDmin);
        trialStiff = grv * pow(dy/dm,alph);
        if(trialStiff < commitStiff)
          trialStiff = commitStiff;
        trialForce = commitForce + trialStiff * (trialDeform - commitDeform);
        trialgss = trialStiff;
        trialStg = 4;
      }

      else {
        if(trialDmin >= commitDeform)
          trialDmin = commitDeform;
        if(trialFmin >= commitForce)
          trialFmin = commitForce;
        dm = fabs(trialDmax);
        if(fabs(trialDmin) > dm)
          dm = fabs(trialDmin);
        trialStiff = grv * pow(dy/dm,alph);
        if(trialStiff < commitStiff)
          trialStiff = commitStiff;
        trialForce = commitForce + trialStiff * (trialDeform - commitDeform);
        trialgss = trialStiff;
        trialStg = 4;
      }
    }

    break;




  //降伏後の初期カーブ上のUmからの除荷
  case 4:

   //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {
     //fumの値の決定
      if(trialForce >  0.0)
        fum = trialFmax;
      if(trialForce <= 0.0)
        fum = trialFmin;

      //再び載荷された状態で除荷点Umを超えない
      if(fabs(trialForce) <= fabs(fum)) {
        trialStiff = commitgss;
        trialForce = commitForce + trialStiff * (trialDeform - commitDeform);
        trialStg = 4;
      }

      //再び載荷された状態で除荷点Umを超える
      else{
        if(trialForce > 0.0) {
          trialStiff = gy;
          trialForce = fy + gy * (trialDeform - dy);
          trialStg = 3;
        }
        else {
          trialStiff = gy;
          trialForce = -fy + gy * (trialDeform + dy);
          trialStg = 3;
        }
      }
    }


    //unloading
    else {


      if(commitForce * trialForce >= 0.0) {

        trialStiff = commitStiff;
        trialForce = trialForce;
        trialStg = 4;

      }
      //load reversal
      else{

        //反対側にクラックが生じていない場合
        if(commitCrk <= 1){

          //荷重がクラック点に達していない
          if(fabs(trialForce) <= fc) {
            trialStiff = commitgss;
            trialForce = commitForce + commitgss * (trialDeform - commitDeform);
            trialStg = 15;
          }

          //荷重がクラック点を超える
          if(fabs(trialForce) > fc) {

            //荷重が降伏点に達していない
            if(fabs(trialForce) <= fy){
              if(trialForce < 0.0) {
                fum = commitFmax;
                dum = commitDmax;
                trialDcNeg = (-fc - fum) / commitgss + dum;
                trialgc0 = (-fy + fc) / (-dy - trialDcNeg);
                trialStiff = trialgc0;
                trialForce = -fc + trialgc0 * (trialDeform - trialDcNeg);
                trialStg = 16;
                trialCrk = commitCrk + 1;
              }
              else {
                fum = commitFmin;
                dum = commitDmin;
                trialDcPos = (fc - fum) / commitgss + dum;
                trialgc0 = (fy - fc) / (dy - trialDcPos);
                trialStiff = trialgc0;
                trialForce = fc + trialgc0 * (trialDeform - trialDcPos);
                trialStg = 16;
                trialCrk = commitCrk + 1;
              }
            }
            //荷重が降伏点を超える
            else {
              if(trialForce > 0.0) {
                trialStiff = gy;
                trialForce = fy + gy * (trialDeform - dy);
                trialStg = 3;
                trialCrk = commitCrk + 1;
              }
              else {
                trialStiff = gy;
                trialForce = -fy + gy * (trialDeform + dy);
                trialStg = 3;
                trialCrk = commitCrk + 1;
              }
            }
          }
          if(trialForce > 0.0) {
            trialDmax = dy;
            trialFmax = fy;
          }
          if(trialForce <= 0.0) {
            trialDmin = -dy;
            trialFmin = -fy;
          }
        }

        //反対側にクラックが生じている場合
        else {
          if (fabs(trialForce) > fy) {
            if (trialForce > 0.0) {
              trialStiff = gy;
              trialForce = fy + gy * (trialDeform - dy);
              trialStg =  3;
            }
            if (trialForce <= 0.0) {
              trialStiff = gy;
              trialForce = -fy + gy *(trialDeform + dy);
              trialStg = 3;
            }
          }
          else {
            if (trialForce > 0.0) {
              fum = trialFmax;
              dum = trialDmax;
              trialx0 = commitDeform - commitForce / commitStiff;
            }

            if (trialForce <= 0.0) {
              fum = trialFmin;
              dum = trialDmin;
              trialx0 = commitDeform - commitForce / commitStiff;
            }   
            
            if (fabs(fum) >= fy) {
              trialStiff = fum / ( dum - trialx0 );
              trialForce = trialStiff * (trialDeform - trialx0);
              trialStg = 6;
            }

            if(fabs(fum) < fy){
              ga = fum / ( dum - trialx0 );
                
              if (trialForce > 0.0)
                gb  =  fy / ( dy - trialx0);
              if (trialForce <= 0.0)
                gb  = -fy / (-dy - trialx0);

              if ( gb >= ga ) {
                trialStiff = gb;
                trialForce = gb * (trialDeform - trialx0);

                if (trialForce > 0.0 && trialFmax < fy) {
                  trialFmax = fy;
                  trialDmax = dy;
                  }
                if ( trialForce < 0.0 && trialFmin > -fy) {
                  trialFmin = -fy;
                  trialDmin = -dy;
                }
                trialStg = 6;
              }

              else {
                trialStiff = ga;
                trialForce = ga * (trialDeform - trialx0);
                trialStg = 6;
              }
            }
          }
        }
      }
    }

    break;




  //降伏前の初期カーブ上のUmからの除荷
  case 5:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {
      //再び載荷してまだ除荷点Umを超えない
      if (trialForce >= 0.0)
        fum = trialFmax;
      if (trialForce < 0.0)
        fum = trialFmin;
      if(fabs(trialForce) <= fabs(fum)) {
        trialStiff = commitgss;
        trialForce = commitForce + commitgss * (trialDeform - commitDeform);
        trialStg = 5;
      }
      //再び載荷してまだ除荷点Umを超えない
      else{
        if (trialForce > 0.0){
          trialStiff = commitgc0;
          trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
          trialStg = 2;
        }
        else{
          trialStiff = commitgc0;
          trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
          trialStg = 2;
        }
      }
    }

    //unloading
    else {
      if(commitForce * trialForce >= 0.0) {

        trialStiff = commitStiff;
	      trialForce = trialForce;
        trialStg = 5;

      }

    //loading reversal
      else{
        if(commitCrk <= 1) {
          trialStiff = commitgss;
          trialForce = commitForce + commitgss * (trialDeform - commitDeform);
          trialStg = 14;
        }
        if(commitCrk > 1){
          if(trialForce > 0.0) {
            fum = trialFmax;
            dum = trialDmax;
            trialx0 = commitDeform - commitForce / commitStiff;
          }
          if(trialForce <= 0.0) {
            fum = trialFmin;
            dum = trialDmin;
            trialx0 = commitDeform - commitForce / commitStiff;
          }
          if(fabs(fum) >= fy) {
            trialStiff = fum / (dum - trialx0);
            trialForce = trialStiff * (trialDeform - trialx0);
            trialStg = 6;
          }
          else{
            ga = fum / ( dum - trialx0 );                
            if (commitForce < 0.0){
              gb  =  fy / ( dy - trialx0);
            }
            if (commitForce >= 0.0){
              gb  =  fy / ( dy + trialx0);
            }
            if ( gb >= ga ) {

              if (trialForce > 0.0 && trialFmax < fy) {
                trialFmax = fy;
                trialDmax = dy;
              }
              if ( trialForce < 0.0 && trialFmin > -fy){
                trialFmin = -fy;
                trialDmin = -dy;
              }
              trialStiff = gb;
              trialForce = gb * (trialDeform - trialx0);
              trialStg = 6;
            }
            else {
              trialStiff = ga;
              trialForce = ga * (trialDeform - trialx0);
              trialStg = 6;
            }
          }
        }
      }
    }

    break;


  //初期カーブ上の反対側の点Umを目指す載荷
  case 6:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {

      if(trialForce >= 0.0) {
        fum = trialFmax;
        dum = trialDmax;
      }
      if(trialForce < 0.0) {
        fum = trialFmin;
        dum = trialDmin;
      }
      //載荷してまだUmを超えない
      if(fabs(trialForce) <= fabs(fum)) {
        trialStiff = fum / (dum - commitx0);
        trialForce = trialStiff * (trialDeform - commitx0);
        trialStg = 6;


      }
      //載荷してUmに達する
      if(fabs(trialForce) > fabs(fum)) {

        if(fabs(fum) < fy) {
          if(trialForce > 0.0) {
            trialStiff = commitgc0;
            trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
            trialStg = 2;
          }
          else {
            trialStiff = commitgc0;
            trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
            trialStg = 2;
          }
        }
        else {
          if(trialForce > 0.0) {
            trialStiff = gy;
            trialForce = fy + gy * (trialDeform - dy);
            trialStg = 3;
          }
          else {
            trialStiff = gy;
            trialForce = -fy + gy * (trialDeform + dy);
            trialStg = 3;
          }
        }
      }
    }
    
    //unloading
    else {

      trialFu0 = commitForce;
      trialDu0 = commitDeform;
      trialStiff = commitgss;
      trialForce = commitForce + commitgss * (trialDeform - commitDeform);
      trialStg = 7;

      //(load reversal)
      if (commitForce * trialForce <= 0.0) {

        trialFu0 = commitForce;
        trialDu0 = commitDeform;

        if(trialForce >= 0.0) {
          dum = trialDmax;
          fum = trialFmax;
        }
        else{
          dum = trialDmin;
          fum = trialFmin;
        }

        trialx1 = commitDeform - commitForce / commitgss;
        trialStiff = fum / (dum - trialx1);
        trialForce = trialStiff * (trialDeform - trialx1);
        trialStg = 8;
      }
    }
    break;


  //case 6の除荷点U0からの載荷
  case 7:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {
      //再び載荷して、まだU0を超えない
      if(fabs(trialForce) <= fabs(trialFu0)) {
        trialStiff = commitgss;
        trialForce = commitForce + commitgss * (trialDeform - commitDeform);
        trialStg = 7;
      }
      else {
        if(trialForce >= 0.0){
          dum = trialDmax;
          fum = trialFmax;
        }
        else {
          dum = trialDmin;
          fum = trialFmin;
        }

        trialStiff = fum / (dum - commitx0);
        trialForce = trialStiff * (trialDeform - commitx0);
        trialStg = 6;
      }
    }
    //unloading
    else {
      trialStiff = commitgss;
      trialForce = commitForce + trialStiff * (trialDeform - commitDeform);
      trialStg = 7;

      if (commitForce * trialForce <= 0.0) {
        //load reversal
        if (trialForce >= 0.0){
          dum = trialDmax;
          fum = trialFmax;
        }
        if (trialForce < 0.0) {
          dum = trialDmin;
          fum = trialFmin;
        }
        //Umにまだ達していない
        if (fabs(trialForce) <= fabs(fum)) {
          trialx1    = commitDeform - commitForce / commitStiff;
          trialStiff = fum / (dum - trialx1);
          trialForce = trialStiff * (trialDeform - trialx1);
          trialStg = 8;
        }
        //Umに達した
        if (fabs(trialForce) > fabs(fum)) {
          if (fabs(fum) < fy) {
            if (trialForce < 0.0) {
              trialStiff = commitgc0;
              trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
              trialStg = 2;
            }
            else {
              trialStiff = commitgc0;
              trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
              trialStg = 2;
            }
          }
          else {
            if (trialForce > 0.0) {
              trialStiff = gy;
              trialForce = fy + gy * (trialDeform - dy);
              trialStg = 3;
            }
            else {
              trialStiff = gy;
              trialForce = -fy + gy * (trialDeform + dy);
              trialStg = 3;
            }
          }
        }
      }
    }
    break;

  //初期カーブ上のUmへ向かう載荷
  case 8:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {

    // opserr << "case8_loading" << endln;

      if (trialForce >= 0.0) {
        dum = trialDmax;
        fum = trialFmax;
      }
      if (trialForce < 0.0) {
        dum = trialDmin;
        fum = trialFmin;
      }
      //Umにまだ達していない
      if (fabs(trialForce) <= fabs(fum)) {

        trialStiff = commitStiff;           
	      trialForce = trialForce;
        trialStg = 8;

      }
      //Umに達した
      else {

    // opserr << "case8_2" << endln;

        if (fabs(fum) < fy) {
          if (trialForce > 0.0) {
            trialStiff = commitgc0;
            trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
            trialStg = 2;
          }
          else {
            trialStiff = commitgc0;
            trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
            trialStg = 2;
          }
        }
        else {
          if (trialForce > 0.0) {
            trialStiff = gy;
            trialForce = fy + gy * (trialDeform - dy);
            trialStg = 3;
          }
          else {
            trialStiff = gy;
            trialForce = -fy + gy * (trialDeform + dy);
            trialStg = 3;
          }
        }
      }
    }

    //unloading
    else {

      trialFu1 = commitForce;
      trialDu1 = commitDeform;
      trialStiff = commitgss;
      trialForce = trialFu1 + commitgss * (trialDeform - trialDu1);
      trialStg = 9;



    //(load reversal)
      if(commitForce * trialForce <= 0.0) {

        trialFu1 = commitForce;
        trialDu1 = commitDeform;
        trialx2 = commitDeform - commitForce / commitgss;
        trialStiff = commitFu0 / (commitDu0 - trialx2);
        trialForce = trialStiff * (trialDeform - trialx2);
        trialStg = 10;
      }
    }
    break;

  //case 8の除荷点U1からの除荷
  case 9:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {

      //u1に達するまで
      if (fabs(trialForce) <= fabs(commitFu1)) {

        trialStiff = commitStiff;
        trialForce = trialForce; 
        trialStg = 9;

      }
      //u1に達した
      else {

        if(trialForce >= 0.0) {
          dum = trialDmax;
          fum = trialFmax;
        }
        if(trialForce < 0.0) {
          dum = trialDmin;
          fum = trialFmin;
        }
        if(fabs(trialForce) <= fabs(fum)) {
          if(trialForce >= 0.0) {
            dum = trialDmax;
            fum = trialFmax;
          }
          if(trialForce < 0.0) {
            dum = trialDmin;
            fum = trialFmin;
          }
          trialStiff = fum / (dum - commitx1);
          trialForce = trialStiff * (trialDeform - commitx1);
          trialStg = 8;
        }
        else {
          if(fabs(fum) < fy) {
            if (trialForce > 0.0) {
              trialStiff = commitgc0;
              trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
              trialStg = 2;
            }
            else {
              trialStiff = commitgc0;
              trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
              trialStg = 2;
              }
            }
          else {
            if (trialForce > 0.0) {
              trialStiff = gy;
              trialForce = fy + gy * (trialDeform - dy);
              trialStg = 3;
            }
            else {
              trialStiff = gy;
              trialForce = -fy + gy * (trialDeform + dy);
              trialStg = 3;
            }
          }
        }
      }
    }

    //unloading
    else {

        trialStiff = commitStiff;
        trialForce = trialForce;
        trialStg = 9;

      //load reversal
      if(commitForce * trialForce <= 0.0) {

        //u0に達していない
        if(fabs(trialForce) <= fabs(commitFu0)) {

          trialx2 = commitDeform - commitForce / commitStiff;
          trialStiff = commitFu0 / (commitDu0 - trialx2);
          trialForce = trialStiff * (trialDeform - trialx2);
          trialStg = 10;
        }
        else {
          if(trialForce >= 0.0) {
            fum = trialFmax;
            dum = trialDmax;
          }
          if(trialForce < 0.0) {
            fum = trialFmin;
            dum = trialDmin;
          }
          //Umに達していない
          if(fabs(trialForce) <= fabs(fum)) {
            trialStiff = fum / (dum - commitx0);
            trialForce = trialStiff * (trialDeform - commitx0);
            trialStg = 6;
          }
          //Umに達した
          else {
            if(fabs(fum) < fy) {
              if (trialForce > 0.0) {
                trialStiff = commitgc0;
                trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
                trialStg = 2;
              }
              else {
                trialStiff = commitgc0;
                trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
                trialStg = 2;
              }
            }
            else {
              if (trialForce > 0.0) {
                trialStiff = gy;
                trialForce = fy + gy * (trialDeform - dy);
                trialStg = 3;
              }
              else {
                trialStiff = gy;
                trialForce = -fy + gy * (trialDeform + dy);
                trialStg = 3;
              }
            }
          }
        }
      }
    }
    break;


  //case 9から引き続き、荷重ゼロを越えて、反対側のU0へ向かう
  case 10:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {

      //u0に達していない
      if(fabs(trialForce) <= fabs(commitFu0)) {

	      trialStiff = commitStiff;
	      trialForce = trialForce;
        trialStg = 10;

      }
      //U0に達した
      else {
        if(trialForce >= 0.0) {
          fum = trialFmax;
          dum = trialDmax;
        }
        if(trialForce < 0.0) {
          fum = trialFmin;
          dum = trialDmin;
        }
        //Umに達していない
        if(fabs(trialForce) <= fabs(fum)) {
          trialStiff = fum / (dum - commitx0);
          trialForce = trialStiff * (trialDeform - commitx0);
          trialStg = 6;
        }
        //Umに達した
        else {
          if(fabs(fum) < fy) {
            if (trialForce > 0.0) {
              trialStiff = commitgc0;
              trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
              trialStg = 2;
            }
            else {
              trialStiff = commitgc0;
              trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
              trialStg = 2;
            }
          }
          else {
            if (trialForce > 0.0) {
              trialStiff = gy;
              trialForce = fy + gy * (trialDeform - dy);
              trialStg = 3;
            }
            else {
              trialStiff = gy;
              trialForce = -fy + gy * (trialDeform + dy);
              trialStg = 3;
            }
          }
        }
      }
	  }
    //unloading
    else{
      trialFu2 = commitForce;
      trialDu2 = commitDeform;
      trialStiff = commitgss;
      trialForce = trialFu2 + commitgss * (trialDeform - trialDu2);
      trialStg = 11;

      if(commitForce * trialForce <= 0.0){

        trialFu2 = commitForce;
        trialDu2 = commitDeform;
        trialx3 = commitDeform - commitForce / commitgss;
        trialStiff = commitFu1 / (commitDu1 - trialx3);
        trialForce = trialStiff * (trialDeform - trialx3);
        trialStg = 12;
      }
    }
    break;

  //case 10の除荷点U2からの除荷
  case 11:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {
      //u2に達していない
      if(fabs(trialForce) <= fabs(commitFu2)) {

        trialStiff = commitStiff;
	      trialForce = trialForce;
        trialStg = 11;

      }
      //u2に達した
      else {

        if(fabs(trialForce) <= fabs(commitFu0)) {
          trialStiff = commitFu0 / (commitDu0 - commitx2);
          trialForce = trialStiff * (trialDeform - commitx2);
          trialStg = 10;
        }
        //U0に達した
        else {
          if(trialForce >= 0.0) {
            fum = trialFmax;
            dum = trialDmax;
          }
          if(trialForce < 0.0) {
            fum = trialFmin;
            dum = trialDmin;
          } 
          //Umに達していない
          if(fabs(trialForce) <= fabs(fum)) {
            trialStiff = fum / (dum - commitx0);
            trialForce = trialStiff * (trialDeform - commitx0);
            trialStg = 6;
          }
          //Umに達した
          else {
            if(fabs(fum) < fy) {
              if (trialForce > 0.0) {
                trialStiff = commitgc0;
                trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
                trialStg = 2;
              }
              else {
                trialStiff = commitgc0;
                trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
                trialStg = 2;
              }
            }
            else {
              if (trialForce > 0.0) {
                trialStiff = gy;
                trialForce = fy + gy * (trialDeform - dy);
                trialStg = 3;
              }
              else {
                trialStiff = gy;
                trialForce = -fy + gy * (trialDeform + dy);
                trialStg = 3;
              }
            }
          }
        }
      }
    }
    //unloading
    else {

      trialStiff = commitStiff;
      trialForce = trialForce;
      trialStg = 11;

      //load reversal
      if(commitForce * trialForce <= 0.0) {

        //u1に達していない
        if(fabs(trialForce) <= fabs(commitFu1)) {

          trialx3 = commitDeform - commitForce / commitStiff;
          trialStiff = commitFu1 / (commitDu1 - trialx3);
          trialForce = trialStiff * (trialDeform - trialx3);
          trialStg = 12;
        }
        //u1に達した
        else {
          if(trialForce >= 0.0) {
            fum = trialFmax;
            dum = trialDmax;
          }
          if(trialForce < 0.0) {
            fum = trialFmin;
            dum = trialDmin;
          }
          //Umに達していない
          if(fabs(trialForce) <= fabs(fum)) {
            trialStiff = fum / (dum - commitx1);
            trialForce = trialStiff * (trialDeform - commitx1);
            trialStg = 8;
          }
          //Umに達した
          else {
            if(fabs(fum) < fy) {
              if (trialForce > 0.0) {
                trialStiff = commitgc0;
                trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
                trialStg = 2;
              }
              else {
                trialStiff = commitgc0;
                trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
                trialStg = 2;
              }
            }
            else {
              if (trialForce > 0.0) {
                trialStiff = gy;
                trialForce = fy + gy * (trialDeform - dy);
                trialStg = 3;
              }
              else {
                trialStiff = gy;
                trialForce = -fy + gy * (trialDeform + dy);
                trialStg = 3;
              }
            }
          }
        }
      }
    }
    break;


  //case 11から引き続き、荷重を越えて、反対側のU1へ向かう
  case 12:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {

      //u1に達していない
      if(fabs(trialForce) <= fabs(commitFu1)) {

        trialStiff = commitStiff;
        trialForce = trialForce;
        trialStg = 12;

      }
      
      //U1に達した
      else {
        if(trialForce >= 0.0) {
          fum = trialFmax;
          dum = trialDmax;
        }
        if(trialForce < 0.0) {
          fum = trialFmin;
          dum = trialDmin;
        } 
        //Umに達していない
        if(fabs(trialForce) <= fabs(fum)) {
          trialStiff = fum / (dum - commitx1);
          trialForce = trialStiff * (trialDeform - commitx1);
          trialStg = 8;
        }


        //Umに達した
        else {
          if(fabs(fum) < fy) {
            if (trialForce > 0.0) {
              trialStiff = commitgc0;
              trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
              trialStg = 2;
            }
            else {
              trialStiff = commitgc0;
              trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
              trialStg = 2;
            }
          }
          else {
            if (trialForce > 0.0) {
              trialStiff = gy;
              trialForce = fy + gy * (trialDeform - dy);
              trialStg = 3;
            }
            else {
              trialStiff = gy;
              trialForce = -fy + gy * (trialDeform + dy);
              trialStg = 3;
            }
          }
        }
      }
    }
    //unloading
    else {
      trialFu3 = commitForce;
      trialDu3 = commitDeform;
      trialStiff = commitgss;
      trialForce = commitForce + commitgss * (trialDeform - trialDu3);
      trialStg = 13;

      //load reversal
      if(commitForce * trialForce <= 0.0) {

        trialFu3 = commitStiff;
        trialDu3 = commitDeform;
        trialx2 = commitDeform - commitForce / commitgss;
        trialStiff = commitFu0 / (commitDu0 - trialx2);
        trialForce = trialStiff * (trialDeform - trialx2);
        trialStg = 10;
      }
    }

    break;

  //case 12の除荷点U3からの除荷
  case 13:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {
      //U3に達していない
      if(fabs(trialForce) <= fabs(commitFu3)){

        trialStiff = commitStiff;
        trialForce = trialForce;
        trialStg = 13;

      }
      //U3に達した
      else {

        //u1に達していない
        if(fabs(trialForce) <= fabs(commitFu1)) {
          trialStiff = commitFu1 / (commitDu1 - commitx3);
          trialForce = trialStiff * (trialDeform - commitx3);
          trialStg = 12;
        }

        //u1に達した
        else {
          if(trialForce >= 0.0) {
            fum = trialFmax;
            dum = trialDmax;
          }
          if(trialForce < 0.0) {
            fum = trialFmin;
            dum = trialDmin;
          }
          //Umに達していない
          if(fabs(trialForce) <= fabs(fum)) {
            trialStiff = fum / (dum - commitx1);
            trialForce = trialStiff * (trialDeform - commitx1);
            trialStg = 8;
          }
          //Umに達した
          else {
            if(fabs(fum) < fy) {
              if (trialForce > 0.0) {
                trialStiff = commitgc0;
                trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
                trialStg = 2;
              }
              else {
                trialStiff = commitgc0;
                trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
                trialStg = 2;
              }
            }
            else {
              if (trialForce > 0.0) {
                trialStiff = gy;
                trialForce = fy + gy * (trialDeform - dy);
                trialStg = 3;
              }
              else {
                trialStiff = gy;
                trialForce = -fy + gy * (trialDeform + dy);
                trialStg = 3;
              }
            }
          }
        }
      }
    }
    //unloading
    else {

        trialStiff = commitStiff;
        trialForce = trialForce;
        trialStg = 13;

      //load reversal
      if(commitForce * trialForce <= 0.0) {

        //u0に達していない
        if(fabs(trialForce) <= fabs(commitFu0)) {
          trialx2 = commitDeform - commitForce / commitStiff;
          trialStiff = commitFu0 / (commitDu0 - trialx2);
          trialForce = trialStiff * (trialDeform - trialx2);
          trialStg = 10;
        }
        //u0に達した
        else {
          if(trialForce >= 0.0) {
            fum = trialFmax;
            dum = trialDmax;
          }
          if(trialForce < 0.0) {
            fum = trialFmin;
            dum = trialDmin;
          }
          //Umに達していない
          if(fabs(trialForce) <= fabs(fum)) {
            trialStiff = fum / (dum - commitx0);
            trialForce = trialStiff * (trialDeform - commitx0);
            trialStg = 6;
          }
          //Umに達した
          else {
            if(fabs(fum) < fy) {
              if (trialForce > 0.0) {
                trialStiff = commitgc0;
                trialForce = fc + commitgc0 * (trialDeform - commitDcPos);
                trialStg = 2;
              }
              else {
                trialStiff = commitgc0;
                trialForce = -fc + commitgc0 * (trialDeform - commitDcNeg);
                trialStg = 2;
              }
            }
            else {
              if (trialForce > 0.0) {
                trialStiff = gy;
                trialForce = fy + gy * (trialDeform - dy);
                trialStg = 3;
              }
              else {
                trialStiff = gy;
                trialForce = -fy + gy * (trialDeform + dy);
                trialStg = 3;
              }
            }
          }
        }
      }
    }
    break;

  //case 5から引き続き、一方でクラックが生じてから、他方のクラックの入っていない方向への載荷
  case 14:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {


      //他方のクラック点を越えるまで
      if(fabs(trialForce) <= fc) {
        trialStiff = commitgss;
        trialForce = commitForce + commitgss * (trialDeform - commitDeform);
        trialStg = 14;
      }



      //他方のクラック点を越えたら
      else {
        if (trialForce > 0.0) {
          trialStiff = commitgc0;
          trialForce = fc + commitgc0 * (trialDeform - dc);
          trialStg = 2;
          trialCrk = commitCrk + 1;
        }
        else {
          trialStiff = commitgc0;
          trialForce = -fc + commitgc0 * (trialDeform + dc);
          trialStg = 2;
          trialCrk = commitCrk + 1;
        }
      }
    }
    //unloading
    else {
      trialStiff = commitgss;
      trialForce = commitForce + commitgss * (trialDeform - commitDeform);
      trialStg = 14;


      //load reversal
      if(commitForce * trialForce <= 0.0) {
        trialStiff = commitgss;
        trialForce = commitForce + commitgss * (trialDeform - commitDeform);
        trialStg = 5;
      }
    }

    break;



  //case 4から引き続き、一方で降伏してから、他方のクラックの入っていない方向への載荷
  case 15:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {

      //他方のクラック点を越えるまで
      if(fabs(trialForce) <= fc) {
        trialStiff = commitgss;
        trialForce = commitForce + commitgss * (trialDeform - commitDeform);
        trialStg = 15;
      }
      //他方のクラック点を越えたら
      else {
        if(trialForce < 0.0) {
          fum = commitFmax;
          dum = commitDmax;
          trialDcNeg = (-fc-fum) / commitgss + dum;
          trialgc0 = (-fy + fc) / (-dy-trialDcNeg);
          trialStiff = trialgc0;
          trialForce = -fc + trialgc0 * (trialDeform - trialDcNeg);
          trialStg = 16;
          trialCrk = commitCrk + 1;
        }
        else {
          fum = commitFmin;
          dum = commitDmin;
          trialDcPos = (fc-fum) / commitgss + dum;
          trialgc0 = (fy-fc) / (dy-trialDcPos);
          trialStiff = trialgc0;
          trialForce = fc + trialgc0 * (trialDeform - trialDcPos);
          trialStg = 16;
          trialCrk = commitCrk + 1;
        }
      }
    }

    //unloading
    else {
      trialStiff = commitgss;
      trialForce = commitForce + commitgss * (trialDeform - commitDeform);
      trialStg = 15;

      //load reversal
      if(commitForce * trialForce <= 0.0) {
        if(trialForce >= 0.0) {
          fum = trialFmax;
          dum = trialDmax;
        }
        if(trialForce < 0.0) {
          fum = trialFmin;
          dum = trialDmin;
        }
        //Umに達していない
        if(fabs(trialForce) <= fabs(fum)) {
          trialStiff = commitgss;
          trialForce = commitForce + commitgss * (trialDeform - commitDeform);
          trialStg = 4;
        }
        //Umに達した
        else {
          if (trialForce > 0.0) {
            trialStiff = gy;
            trialForce = fy + gy * (trialDeform - dy);
            trialStg = 3;
          }
          else {
            trialStiff = gy;
            trialForce = -fy + gy * (trialDeform + dy);
            trialStg = 3;
          }
        }
      }
    }

    break;



  //case 15から引き続き、降伏点をめざす載荷
  case 16:
    //loading
    if((commitForce >= 0.0 && dDeform >= 0.0) || (commitForce <= 0.0 && dDeform <= 0.0)) {
      //降伏点を越えるまで
      if(fabs(trialForce) <= fy) {

        trialStiff = commitStiff;
        trialForce = trialForce;
	      trialStg = 16;
      }
      else {
        if (trialForce > 0.0) {
          trialStiff = gy;
          trialForce = fy + gy * (trialDeform - dy);
          trialStg = 3;
        }
        else {
          trialStiff = gy;
          trialForce = -fy + gy * (trialDeform + dy);
          trialStg = 3;
        }
      }
    }
    //unloading
    else {
      trialFu0 = commitForce;
      trialDu0 = commitDeform;
      trialx0 = trialDu0 - trialFu0 / commitStiff;
      trialStiff = commitgss;
      trialForce = commitForce + commitgss * (trialDeform - trialDu0);
      trialStg = 7;

      //load reversal
      if(commitForce * trialForce <= 0.0) {
        trialFu0 = commitForce;
        trialDu0 = commitDeform;
        trialx0 = trialDu0 - trialFu0 / commitStiff;
        trialx1 = trialDu0 - trialFu0 / commitgss;

        if(trialForce > 0.0) {
          fum = trialFmin;
          dum = trialDmax;
        }
        else {
          fum = trialFmin;
          dum = trialDmin;
        }
        trialStiff = fum / (dum - trialx1);
        trialForce = trialStiff * (trialDeform - trialx1);
        trialStg = 8;
      }
    }

    break;
  }

  // check print--------------------------------------------
   //opserr << "---Takeda::setTrialStrain end\n";
   //opserr << " trialStg = " << trialStg <<endln;
   //opserr << " commitStg = " << commitStg <<endln;
  // -------------------------------------------------------

return 0;

}



double 
Takeda::getStress(void)
{
  return trialForce;
}

double 
Takeda::getTangent(void)
{
  return trialStiff;
}

double 
Takeda::getInitialTangent(void)
{
  return ge;
}


double 
Takeda::getStrain(void)
{
  return trialDeform;
}

int 
Takeda::commitState(void)
{


  // check print--------------------------------------------
  // opserr << "---Takeda::commitStateA\n";
  // opserr << "---commitState2_commitStg = " << commitStg <<endln;
  // -------------------------------------------------------

  commitDeform = trialDeform;
  commitForce  = trialForce;
  commitStiff  = trialStiff;

  commitStg   =  trialStg;
  commitCrk   =  trialCrk;
  commitDmax  =  trialDmax;
  commitFmax  =  trialFmax;
  commitDmin  =  trialDmin;
  commitFmin  =  trialFmin;
  commitDu0   =  trialDu0;
  commitFu0   =  trialFu0;
  commitDu1   =  trialDu1;
  commitFu1   =  trialFu1;
  commitDu2   =  trialDu2;
  commitFu2   =  trialFu2;
  commitDu3   =  trialDu3;
  commitFu3   =  trialFu3;
  commitx0    =  trialx0;
  commitx1    =  trialx1;
  commitx2    =  trialx2;
  commitx3    =  trialx3;
  commitgss   =  trialgss;
  commitDcPos =  trialDcPos;
  commitDcNeg =  trialDcNeg;
  commitgc0   =  trialgc0;


  // opserr << trialStg << endln;
// opserr << trialgc0 << endln;
  // check print--------------------------------------------
  // opserr << "---Takeda::commitStateB\n";
  // opserr << "---commitState2_commitStg = " << commitStg <<endln;
  // -------------------------------------------------------
  return 0;
}


int 
Takeda::revertToLastCommit(void)
{

  trialDeform  = commitDeform;
  trialForce   = commitForce;
  trialStiff   = commitStiff;

  trialStg    =  commitStg;
  trialCrk    =  commitCrk;
  trialDmax   =  commitDmax;
  trialFmax   =  commitFmax;
  trialDmin   =  commitDmin;
  trialFmin   =  commitFmin;
  trialDu0    =  commitDu0;
  trialFu0    =  commitFu0;
  trialDu1    =  commitDu1;
  trialFu1    =  commitFu1;
  trialDu2    =  commitDu2;
  trialFu2    =  commitFu2;
  trialDu3    =  commitDu3;
  trialFu3    =  commitFu3;
  trialx0     =  commitx0;
  trialx1     =  commitx1;
  trialx2     =  commitx2;
  trialx3     =  commitx3;
  trialgss    =  commitgss;
  trialDcPos  =  commitDcPos;
  trialDcNeg  =  commitDcNeg;
  trialgc0    =  commitgc0;


  // check print--------------------------------------------
  // opserr << "---Takeda::revertToLastCommit\n";
  // -------------------------------------------------------

  return 0;
}

int 
Takeda::revertToStart(void)
{

  trialDeform  = 0.0;
  trialForce   = 0.0;
  trialStiff   = ge;
  commitDeform = 0.0;
  commitForce  = 0.0;
  commitStiff  = ge;


  dm   = 0.0;
  ga   = 0.0;
  gb   = 0.0;
  dum  = 0.0;
  fum  = 0.0;
  gc   = bc * ge;
  gy   = by * ge;
  dc   = fc / ge;
  dy   = dc + (fy - fc) / gc;
  grv = (fc + fy) / (dc + dy);


  trialStg    = 1;
  trialCrk    = 0;
  trialDmax   = 0.0;
  trialFmax   = 0.0;
  trialDmin   = 0.0;
  trialFmin   = 0.0;

  trialDu0    = 0.0;
  trialFu0    = 0.0;
  trialDu1    = 0.0;
  trialFu1    = 0.0;
  trialDu2    = 0.0;
  trialFu2    = 0.0;
  trialDu3    = 0.0;
  trialFu3    = 0.0;

  trialx0     = 0.0;
  trialx1     = 0.0;
  trialx2     = 0.0;
  trialx3     = 0.0;
  trialgss    = 0.0;
  trialDcPos  =  dc;
  trialDcNeg  = -dc;
  trialgc0    =  gc;

  commitStg   = 1;
  commitCrk   = 0;
  commitDmax  = 0.0;
  commitFmax  = 0.0;
  commitDmin  = 0.0;
  commitFmin  = 0.0;
  commitDu0   = 0.0;
  commitFu0   = 0.0;
  commitDu1   = 0.0;
  commitFu1   = 0.0;
  commitDu2   = 0.0;
  commitFu2   = 0.0;
  commitDu3   = 0.0;
  commitFu3   = 0.0;
  commitx0    = 0.0;
  commitx1    = 0.0;
  commitx2    = 0.0;
  commitx3    = 0.0;
  commitgss   = 0.0;
  commitDcPos =  dc;
  commitDcNeg = -dc;
  commitgc0   =  gc;

  // check print--------------------------------------------
  //  opserr << "---Takeda::revertToStart\n";
  // -------------------------------------------------------

  return 0;
}


UniaxialMaterial *
Takeda::getCopy(void)
{
  Takeda *theCopy = new Takeda(this->getTag(), fc, fy, ge, 
					 bc, by, alph);

  // Copy 
  theCopy->trialDeform  = trialDeform;
  theCopy->trialForce   = trialForce;
  theCopy->trialStiff   = trialStiff;
  theCopy->commitDeform = commitDeform;
  theCopy->commitForce  = commitForce;
  theCopy->commitStiff  = commitStiff;

  theCopy->dm   = dm;
  theCopy->ga   = ga;
  theCopy->gb   = gb;
  theCopy->dum  = dum;
  theCopy->fum  = fum;
  theCopy->gc   = gc;
  theCopy->gy   = gy;
  theCopy->dc   = dc;
  theCopy->dy   = dy;
  theCopy->grv  = grv;

  theCopy->trialStg    = trialStg;
  theCopy->trialCrk    = trialCrk;
  theCopy->trialDmax   = trialDmax;
  theCopy->trialFmax   = trialFmax;
  theCopy->trialDmin   = trialDmin;
  theCopy->trialFmin   = trialFmin;

  theCopy->trialDu0    = trialDu0;
  theCopy->trialFu0    = trialFu0;
  theCopy->trialDu1    = trialDu1;
  theCopy->trialFu1    = trialFu1;
  theCopy->trialDu2    = trialDu2;
  theCopy->trialFu2    = trialFu2;
  theCopy->trialDu3    = trialDu3;
  theCopy->trialFu3    = trialFu3;

  theCopy->trialx0     = trialx0;
  theCopy->trialx1     = trialx1;
  theCopy->trialx2     = trialx2;
  theCopy->trialx3     = trialx3;
  theCopy->trialgss    = trialgss;
  theCopy->trialDcPos  = trialDcPos;
  theCopy->trialDcNeg  = trialDcNeg;
  theCopy->trialgc0    = trialgc0;

  theCopy->commitStg   = commitStg;
  theCopy->commitCrk   = commitCrk;
  theCopy->commitDmax  = commitDmax;
  theCopy->commitFmax  = commitFmax;
  theCopy->commitDmin  = commitDmin;
  theCopy->commitFmin  = commitFmin;
  theCopy->commitDu0   = commitDu0;
  theCopy->commitFu0   = commitFu0;
  theCopy->commitDu1   = commitDu1;
  theCopy->commitFu1   = commitFu1;
  theCopy->commitDu2   = commitDu2;
  theCopy->commitFu2   = commitFu2;
  theCopy->commitDu3   = commitDu3;
  theCopy->commitFu3   = commitFu3;
  theCopy->commitx0    = commitx0;
  theCopy->commitx1    = commitx1;
  theCopy->commitx2    = commitx2;
  theCopy->commitx3    = commitx3;
  theCopy->commitgss   = commitgss;
  theCopy->commitDcPos = commitDcPos;
  theCopy->commitDcNeg = commitDcNeg;
  theCopy->commitgc0   = commitgc0;

  // check print--------------------------------------------
  //  opserr << "---Takeda::getCopy\n";
  // -------------------------------------------------------
  return theCopy;
}


int 
Takeda::sendSelf(int cTag, Channel &theChannel)
{

  int res = 0;


  static Vector data(67);

  data(0)  = this->getTag();

  data(1)  = fc;
  data(2)  = fy;
  data(3)  = ge;
  data(4)  = bc;
  data(5)  = by;
  data(6)  = alph;

  data(7)  = trialDeform;
  data(8)  = trialForce;
  data(9)  = trialStiff;
  data(10)  = commitDeform;
  data(11)  = commitForce;
  data(12)  = commitStiff;

  data(13)  = dm;
  data(14)  = ga;
  data(15)  = gb;
  data(16)  = dum;
  data(17)  = fum;
  data(18)  = gc;
  data(19)  = gy;
  data(20)  = dc;
  data(21)  = dy;
  data(22)  = grv;

  data(23)  = trialStg;
  data(24)  = trialCrk;
  data(25)  = trialDmax;
  data(26)  = trialFmax;
  data(27)  = trialDmin;
  data(28)  = trialFmin;
  data(29)  = trialDu0;
  data(30)  = trialFu0;
  data(31)  = trialDu1;
  data(32)  = trialFu1;
  data(33)  = trialDu2;
  data(34)  = trialFu2;
  data(35)  = trialDu3;
  data(36)  = trialFu3;
  data(37)  = trialx0;
  data(38)  = trialx1;
  data(39)  = trialx2;
  data(40)  = trialx3;
  data(41)  = trialgss;
  data(42)  = trialDcPos;
  data(43)  = trialDcNeg;
  data(44)  = trialgc0;

  data(45)  = commitStg;
  data(46)  = commitCrk;
  data(47)  = commitDmax;
  data(48)  = commitFmax;
  data(49)  = commitDmin;
  data(50)  = commitFmin;
  data(51)  = commitDu0;
  data(52)  = commitFu0;
  data(53)  = commitDu1;
  data(54)  = commitFu1;
  data(55)  = commitDu2;
  data(56)  = commitFu2;
  data(57)  = commitDu3;
  data(58)  = commitFu3;
  data(59)  = commitx0;
  data(60)  = commitx1;
  data(61)  = commitx2;
  data(62)  = commitx3;
  data(63)  = commitgss;
  data(64)  = commitDcPos;
  data(65)  = commitDcNeg;
  data(66)  = commitgc0;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Takeda::sendSelf() - failed to send data\n";

  // check print--------------------------------------------
  //  opserr << "---Takeda::sendSelf\n";
  // -------------------------------------------------------
  return res;
}

int 
Takeda::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(67);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);

  if (res < 0) {
      opserr << "Takeda::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }


  else {
    this->setTag((int)data(0));

  fc            = data(1);
  fy            = data(2);
  ge            = data(3);
  bc            = data(4);
  by            = data(5);
  alph          = data(6);

  trialDeform   = data(7);
  trialForce    = data(8);
  trialStiff    = data(9);
  commitDeform  = data(10);
  commitForce   = data(11);
  commitStiff   = data(12);

  dm            = data(13);
  ga            = data(14);
  gb            = data(15);
  dum           = data(16);
  fum           = data(17);
  gc            = data(18);
  gy            = data(19);
  dc            = data(20);
  dy            = data(21);
  grv          = data(22);


  trialStg    = (int)data(23);
  trialCrk    = (int)data(24);
  trialDmax   = data(25);
  trialFmax   = data(26);
  trialDmin   = data(27);
  trialFmin   = data(28);
  trialDu0    = data(29);
  trialFu0    = data(30);
  trialDu1    = data(31);
  trialFu1    = data(32);
  trialDu2    = data(33);
  trialFu2    = data(34);
  trialDu3    = data(35);
  trialFu3    = data(36);
  trialx0     = data(37);
  trialx1     = data(38);
  trialx2     = data(39);
  trialx3     = data(40);
  trialgss    = data(41);
  trialDcPos  = data(42);
  trialDcNeg  = data(43);
  trialgc0    = data(44);

  commitStg   = (int)data(45);
  commitCrk   = (int)data(46);
  commitDmax  = data(47);
  commitFmax  = data(48);
  commitDmin  = data(49);
  commitFmin  = data(50);
  commitDu0   = data(51);
  commitFu0   = data(52);
  commitDu1   = data(53);
  commitFu1   = data(54);
  commitDu2   = data(55);
  commitFu2   = data(56);
  commitDu3   = data(57);
  commitFu3   = data(58);
  commitx0    = data(59);
  commitx1    = data(60);
  commitx2    = data(61);
  commitx3    = data(62);
  commitgss   = data(63);
  commitDcPos = data(64);
  commitDcNeg = data(65);
  commitgc0   = data(66);
  }


  // check print--------------------------------------------
  // opserr << "---Takeda::recvSelf\n";
  // -------------------------------------------------------
  return res;
}

void 
Takeda::Print(OPS_Stream &s, int flag)
{
  s << "Takeda : " << this->getTag() << endln;
  s << "  fc: " << fc << endln;
  s << "  fy: " << fy << endln;
  s << "  ge: " << ge << endln;
  s << "  bc: " << bc << endln;
  s << "  alph: " << alph << endln;
  s << "  by: " << by << endln;

  return;
}
