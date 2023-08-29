#include <mathsym.h>
///\file globalvariables.cc
///\brief Global variables
///
///Global variables used in LignumForest initialized by GrowthLoop.
///Variables are dummy in the context of CrownDensity, needed for compilation.
namespace LignumForest{
  ///\defgroup spaceb Global variables to set heights and lengths
  ///@{
  ///Initial height ot trees \sa pine-em98.L
  double H_0_ini = 0.0;
  ///Variation of initial tree heights
  double H_var_ini = 0.0;
  
  ///Variation in the initial max number of buds.
  int n_buds_ini_max = 0.0;
  ///Variation in the initial min number of buds.
  int n_buds_ini_min = 0.0; 
  ///If bud variation is on \sa rel_bud
  bool bud_variation = false;
  ///The effect of crowding on number of lateral buds via function bud_view_f.
  ///Requires \sa bud_variation
  double rel_bud = 0.0;
  ///\deprecated Variation in branch angle. Set in \sa pine-em98.L
  double branch_angle = 45.0 * PI_VALUE / 180.0;
}
