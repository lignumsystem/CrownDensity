#include <CrownDensityGlobals.h>

namespace Pine{
  double L_age = 0.0;
  double L_H = 0.0;
  const double l1 = 0.9;
}

namespace CrownDensity{  
  ///Random number generator e.g. in L-system
  int ran3_seed=-1234567;
  ///This is needed in SomeFunctors.h SetScotsPineSegmentApical
  bool is_by_branches = false;          

  ///\defgroup adhoc Ad Hoc experiment for segment growth
  ///@{
  ///Increase growth in lower parts of crown
  ParametricCurve adhoc("adhoc.fun");   
  bool is_adhoc = false;
  ParametricCurve fip_mode("fip1.fun");
  ParametricCurve fgo_mode("fgo1.fun");
  ///Random component in segment length in branches.
  ///If true then the effect is on from the throughout the simulation  
  bool is_random_variation=false;
  ///@}

  ///This global variable is needed e.g. in bybranches.cc
  int growthloop_ebh_mode = 0;          

  ///Radiation use efficiency experiment
  LGMdouble max_rueqin;

  ///\defgroup modechange Growth mode change experiment
  ///Architectural and morphological modes now separated and set independently.
  ///`is_mode_change` sets the morphological mode change (reset parameters, functions)
  ///and `is_architecture_change` sets the architecture mode change (L-system implementation) 
  ///@{
  ///This global variable is for change in morphological development
  bool is_mode_change = false;
  ///This global variable is to set year for change in morphological development.
  ///INT_MAX means default is off and must be explicitly set in command line.
  int mode_change_year = INT_MAX;
  ///This global variable denotes change in architecture development
  bool is_architecture_change = false;
  ///This global variable denotes year for change in architecture development
  ///MAX_INT default means off and must be explicitly set in command line.
  int architecture_change_year = INT_MAX;
  ///@}
  ///\defgroup vsoccupance VoxelSpace for space occupancy
  ///@{
  Firmament dummy_firm;
  VoxelSpace space_occupancy(Point(0.0,0.0,0.0),Point(1.0,1.0,1.0),
			     0.1,0.1,0.1,5,5,5,dummy_firm);
  ///@}
  ///\file crowndensity_globals.cc
  ///\brief Globals variables declared and used in CrownDensity

  ///\defgroup hcbconvey For conveying height of grown base to SetScotsPineSegmentLength
  ///@{
  double global_hcb;      
  double dDb;
  ///@}

  ///\defgroup spacecolon Space colonialization experiment
  ///Options for SetScotsPineSegmentLength (in ScotsPine.h)
  ///and other related global variables
  ///@{
  bool space0 = false;
  bool space1 = false;
  bool space2 = false;
  double space2_distance = 0.3;
  bool is_forced_height = false;
  bool is_height_function = false;
  ///@}

  ///\defgroup toptax Top tax experiment
  ///@{
  ParametricCurve toptax("toptax.fun"); ///<function to adjust resource distn among branches
  double tax_share = 0.3;        
  ///@}

  ///This global variable conveys the Bud View Function to Lsystem
  /// and if it is in use
  ///\defgroup budview Bud view experiment
  ///@{
  ParametricCurve bud_view_f("bvf.fun");
  bool is_bud_view_function = false;   
  ///@}
}//end namespace
