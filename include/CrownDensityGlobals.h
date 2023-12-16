#ifndef CROWNDENSITYGLOBALS_H
#define CROWNDENSITYGLOBALS_H
///\file CrownDensityGlobals.h
///\brief Global varables in CrownDensity
///\sa crowndensity_globals.cc

#include <LGMSymbols.h>
#include <Point.h>
#include <ParametricCurve.h>
#include <Firmament.h>
#include <VoxelSpace.h>
///\defgroup lsysglobal L-system global variables
///@{
///These global variables are needed in pine-em98.L 
namespace Pine{
  ///Age
  extern double L_age;
  ///Tree height
  extern double L_H;
  ///A guess for segment shortening
  extern const double l1;
  extern ParametricCurve fnbuds;
  extern ParametricCurve fnbudslight;
  extern ParametricCurve fipbud;
  extern bool is_architecture_change;
  extern int architecture_change_year;
}
///@}

//\defgroup lgmforest Global variables used in LignumForest
///@{
///Global variables used in LignumForest initialized by GrowthLoop.
///Variables are dummy in the context of CrownDensity, needed for compilation.
namespace LignumForest{
  extern double H_0_ini;
  extern double H_var_ini;
  extern bool bud_variation;
  extern int n_buds_ini_max;
  extern int n_buds_ini_min;
  extern double rel_bud;
  extern double branch_angle;
}
///@}

namespace CrownDensity{

  extern int ran3_seed;

  extern ParametricCurve adhoc;
  extern bool is_adhoc;

  extern int growthloop_ebh_mode;       

  extern LGMdouble max_rueqin;

  extern bool is_mode_change;
  extern int mode_change_year;    

  extern Firmament dummy_firm;
  extern VoxelSpace space_occupancy;

  //extern double L_age, L_H;

  extern double global_hcb;      
  extern double dDb;

  extern bool space0;
  extern bool space1;
  extern bool space2;
  extern double space2_distance;
  extern bool is_height_function;

  extern ParametricCurve bud_view_f;
  extern bool is_bud_view_function;   
}//end namespace
#endif
