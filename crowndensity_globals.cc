#include <CrownDensityGlobals.h>
///Random number generator e.g. in L-system
int ran3_seed;
///This is needed in SomeFunctors.h SetScotsPineSegmentApical
bool is_by_branches = false;          

///\defgroup adhoc Ad Hoc experiment
///@{
///Increase growth in lower parts of crown
ParametricCurve adhoc("adhoc.fun");   
bool is_adhoc = false;
///@}

///This global variable is needed e.g. in bybranches.cc
int growthloop_ebh_mode = 0;          

///Radiation use efficiency experiment
LGMdouble max_rueqin;

///\defgroup modechange Growth mode change experiment
///@{
///This global variable is for change in morphological development
bool is_mode_change = false;
///This global variable is to set year for change in morphological development
int mode_change_year = 20;    
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
ParametricCurve bud_view_f;
bool is_bud_view_function = false;   
///@}
