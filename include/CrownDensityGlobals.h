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
///and in its derivatives like pine-em98-branch-C.L
namespace Pine{
  ///\brief Age
  extern double L_age;
  ///\brief Tree height
  extern double L_H;
  ///\brief A guess for segment shortening
  extern const double l1;
  ///\brief Number of buds as a function of foliage mass
  extern ParametricCurve fnbuds;
  extern ParametricCurve fnbudslight;
  ///\brief Number of buds as function of relative light
  extern ParametricCurve fipbud;
  ///\brief Set architecture change on/off
  extern bool is_architecture_change;
  ///\brief Set architecture change year
  extern int architecture_change_year;
}
///@}

///\defgroup lgmforest Global variables used in LignumForest
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
///\defgroup crowndensityglobals Global variables in CrownDensity
///@{
///Global variables used in CrownDensity initialized in main.cc
namespace CrownDensity{

  ///\brief Random seed
  extern int ran3_seed;

  ///\brief Random component in segment length (use to apply conditionally any function implementing segment elongation)
  extern bool is_random_length;
  ///\brief Ad hoc lengthening of shoots at crown base. \sa CrownDensity::adoc
  extern bool is_adhoc;
  ///\brief Function for ad hoc lengthening of shoots at crown base. \sa CrownDensity::is_adhoc
  extern ParametricCurve adhoc;
  ///\brief Fip function after mode change
  extern ParametricCurve fip_mode;
  ///\brief Gravelius order function after mode change
  extern ParametricCurve fgo_mode;

  extern int growthloop_ebh_mode;       

  extern LGMdouble max_rueqin;

  ///\brief Set mode change on/off
  extern bool is_mode_change;
  ///\brief Set mode change year
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
///@}
#endif
