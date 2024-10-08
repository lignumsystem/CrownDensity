#include <CrownDensityGlobals.h>

namespace Pine{
  double L_age = 0.0;
  double L_H = 0.0;
  const double l1 = 0.9;
  ParametricCurve fnbuds;
  ParametricCurve fnbudslight;
  ParametricCurve fipbud;
  bool is_architecture_change = false;
  int architecture_change_year = INT_MAX;
}

namespace CrownDensity{  
  ///Random number generator e.g. in L-system
  int ran3_seed=-1234567;

  ///\defgroup adhoc Ad Hoc experiment for segment growth
  ///@{
  ///Increase growth in lower parts of crown
  ParametricCurve adhoc("adhoc.fun");   
  bool is_adhoc = false;
   ///@}

  ///Explicitely use random component in segment length.
  bool is_random_length = false;
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
   ///This global variable denotes year for change in architecture development
  ///MAX_INT default means off and must be explicitly set in command line.
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
  bool is_height_function = false;
  ///@}

  ///This global variable conveys the Bud View Function to Lsystem
  /// and if it is in use
  ///\defgroup budview Bud view experiment
  ///@{
  ParametricCurve bud_view_f("bvf.fun");
  bool is_bud_view_function = false;
  ///@}

  
  int a_change_start = 0;
  double  max_turn_in_architecture_change = 80.0*PI_VALUE/180.0;
  deque<string> metafile_q;
  deque<int> modechange_year_q;
  

 
}//end namespace
