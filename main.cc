///\file main.cc
///\brief Main growth loop.
///
///Main program growth loop generates now only HDF5 files for result analysis.<br>
///Growth loop steps
///+ Initialize global variables
///+ Read command line
///+ Initialize tree, L-system and global variables from command line
///+ Run the simulation
///+ Save simulation data to HDF5 files.
///\sa CrownDensity
///\page cmakefile Compile CrownDensity
///CMakeLists.txt to compile CrownDensity. There are several L-system files to experiment
///with tree architecture. Check LSYSTEMFILE and LSYSTEMSRC variables in CMakeLists.txt before compilation.
///\include CMakeLists.txt
///\page runscript Run CrownDensity with EBH
///Use the following `run-crowndens.sh` script to run `crowndens` with EBH, ad_hoc etc. experiments.
///Create new scripts or edit `crowndens` command line according to simulation needs.
///\include run-crowndens.sh
///\page runscript_arch Run CrownDensity Basic Model
///Use the following `run-crowndens-basic-model.sh` script to run `crowndens`
///with simple model for segment elongation: \f$L=\lambda \times f_{vi} \times f_{go} \times f_{ip}\f$.
///The command line does not use flags for EBH, ad_hoc etc. experiments.
///
///\include run-crowndens-basic-model.sh
///\page lsystem L-system file
///The following L-system file defines architecture for Scots pine.
///Note that now command line can define architecure change when branches
///will always create two subbranches unless growth conditions set the segment length
///to zero (see GROWTH ARCHITECTURE CHANGE comment)
///\include pine-em98.L
///\page lsystemA L-system experiment A
///The following L-system file experiments with tree architecture.
///Introduce concept *physiological age* for buds. *Each branch bud starts with physiological
///age 1 and each time step increase the physiological age by 1*. When the physiological age
///reaches  CrownDensity::architecture_change_year the terminating bud genererates *always* two side branches.
///In this case branching stops if and only if tree parameters and tree functions determine so.
///\sa \ref lsystemB "L-system experiment B" \ref lsystemC "L-system experiment C" 
///
///\include pine-em98-branch-A.L
///\page lsystemB L-system experiment B
///The following L-system file experiments with tree architecture. Based on pine-em98-branch-A.L.
///Introduce concept *physiological age* for buds. *Each branch bud inherits the physiological
///age of the mother bud  and each time step increase the physiological age by 1*. When the physiological age
///reaches  CrownDensity::architecture_change_year the terminating bud genererates *always* two side branches.
///Only Turn, no Roll along mother axis. In this case branching stops if and only if tree parameters
///and tree functions determine so. \sa \ref lsystemA "L-system experiment A"  \ref lsystemC "L-system experiment C"
///
///\include pine-em98-branch-B.L
///\page lsystemC L-system experiment C
///The following L-system file experiments with tree architecture. Based on  pine-em98-branch-B.L.
///Introduce concept *physiological age* for buds. Each branch bud inherits the physiological
///age of the mother bud  and each time step increase the physiological age by 1. When the physiological age
///reaches  CrownDensity::architecture_change_year the terminating bud genererates *always* two side branches.
///Only Turn, no Roll along mother axis. *These side branches have high branching angle*. See `Pine::turn_branch_max` value in pine-em98-branch-C.L.
///In this case branching stops if and only if tree parameters and tree functions determine so.
///\sa \ref lsystemA "L-system experiment A"  \ref lsystemB "L-system experiment B" 
///
///\include pine-em98-branch-C.L
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <Lignum.h>
#include <Bisection.h>
//XML file 
#include <XMLTree.h>
//Include the implementation of the tree segment and bud
#include <ScotsPine.h>
#include <CrownDensityGlobals.h>
//From ../LignumForest/
#include <GrowthLoop.h> 
#include <Shading.h>
//From ../LignumForest project
#include <TreeDataAfterGrowth.h>
#include <CreateHDF5Files.h>
#include <SomeFunctors.h>         //From ../LignumForest/include
#include <DiameterGrowth.h>       //From ../LignumForest/include
#include <RadiationCrownDens.h>   //From ../LignumForest/include
#include <Palubicki_functors.h>   //From ../LignumForest/include
#include <Space.h>                //From ../LignumForest/include

#if defined (__APPLE__) || defined(__MACOSX__)
#include <VisualFunctor.h>
//Impelements VisualizeLGMTree
#include <GLSettings.h>
#include <OpenGLUnix.h>
#include <LGMVisualization.h>
#endif
//Includes all kinds of stuff, turtle graphics etc.
#include <lengine.h>

//and for pine, see also pine9bp.L in lsys.
///\brief L-system in Pine namespace
namespace Pine{
#include <LSystem.h>

}
using namespace LignumForest;

///\brief Forest simulation with single tree.
///
///Single tree growing in a forest plot. The forest plot consists of copies of the tree
///creating shadowing forest in a voxel space. Outside the voxel space
///forest canopy is assumed to be homogenous with the heigth of tree crown length
///at any time of simulation.
namespace CrownDensity{
  ///\par Main program
  ///Note two command line options:
  ///+ -modeChange that reassigns new functions for growth
  ///+ -architectureChange that changes branching behavious implemented in L-system
  ///\snippet{lineno} main.cc Usagex
  ///\internal
  // [Usagex]
  void Usage()
  {
    cout << "Usage:  ./crowndens  <iter>  <metafile>" <<endl;
    cout << "-numParts <parts> -hw <hw_start>" <<endl;
    cout << "[-toFile <Filename>] [-xml <filename>]" <<endl;
    cout << "[-fipdistrib <filename>] [-writeInterval interval]" << endl;
    cout << "-[-seed <num>] [-selfThinning] [-increaseXi <year>]" <<endl;
    cout << "[-treeFile <filename>] [-voxelCalculation <year>]" << endl;
    cout << "[-space0] [-space1] [-space2] [-standFromFile] [-adHoc] [-randomLength]" << endl;
    cout << "[-budViewFunction] [-EBH] -EBH1 <value>]" << endl;
    cout << "[-EBHREDUCTION <value>] [-EBHFINAL <value>] [-EBHInput <int>] [-RUE <value>]" << endl;
    cout << "[-space2Distance <Value>] [-minLeaderLen <value>]" << endl;
    cout << "[-modeChange <year> ] [-architectureChange <year> [-aChangeStart <year>]]" << endl;
    cout << "[-kBorderConifer]" << endl;

    cout << endl;
    cout << "-adHoc             Increase growth in lower parts of the crown." <<endl;
    cout << "EBH resource distn can be in use in two ways. Both are set by command line arguments." << endl;
    cout << "-EBH               means EBH is in use and values (of lambda parameter) are" << endl;
    cout << "                   specified for all Gravelius orders (function SPEBHF, ScotsPine.h, function" << endl;
    cout << "                   file is specified the constructor of the tree)." << endl;
    cout << "-EBH1 <value>      means that EBH is in use and one value <value> is used" << endl;
    cout << "                   for all Gravelius orders. Option -EBH1 <value> overrides option -EBH" << endl;
    cout << "                   (EBH is set as SPis_EBH (Scots Pine Parameter Double SPPD), thus 0 == false, 1 == true)" << endl;
    cout << "                   EBH is according to W. Palubicki and K. Horel and S. Longay and" << endl;
    cout << "                   A. Runions and B. Lane and R. Mech and P. Prusinkiewicz. 2009." << endl;
    cout << "                   Self-organizing tree models for image synthesis ACM Transactions on Graphics 28 58:1-10." << endl;
    cout << "-EBHREDUCTION <value> If values of EBH parameters for all orders are reduced or increased  as" << endl;
    cout << "                      new_value = <value> * prev_value in each year after year 20. The goal of change" << endl;
    cout << "                   can be set by -EBHFINAL. The default value is 0.5"<< endl;
    cout << "-EBHFINAL <value>  Sets the goal value of -EBHREDUCTION" << endl;
    cout << "-EBHInput <int>    Changes the variable that runs EBH. If int == 1, it is Qabs, if int == 2, it is rue*Qabs," << endl;
    cout << "                   any other value or missing -EBHInput means Qin runs EBH." << endl;
    cout << "-RUE <value>       The radiation use effeciency (rue) varies as a function of TreeSegments initial radiation" << endl;
    cout << "                   conditions. Photosynthetic production of TreeSegment = rue * LGApr * Qabs. <value> = degree of" << endl;
    cout << "                   increase of rue as a function of shadiness (0 < <value> < 2)." << endl;
    cout << "-modeChange <year> If the mode of morphological development changes in year <year>" << endl;
    cout << "                   In this case morphology changes to fip.fun and fgo.fun (may be e.g. EBH before)" << endl;
    cout << "                   and new functions for fip (fip1.fun) and fgo (fgo1.fun) are read in." << endl;
    cout << "-randomLength      Use explicitely random component in segment length." << endl;
    cout << "-architectureChange <year> If the mode of architectural development in a shoot (Segment) changes" << endl;
    cout << "                   when it is <year> years old. This is implemented in the L system." <<endl;
    cout << "-aChangeStart <year>    If -architectureChange the time when -architectureChange starts to work. Default" << endl;
    cout << "                    value is 0. Note that -architectureChange will" << endl;
    cout << "                    affect morphological development when a shoot (Segment) is <year> years old. In order" << endl;
    cout << "                    to start the change in morphology in year tm, <year> in -aChangeStart must be set as" << endl;
    cout << "                    tm - <year> of -architectureChange." << endl;
    cout << "                    See the global variables is_architecure_change and architecture_change_year" <<endl;       
    cout << "-kBorderConifer <value>   Extinction coefficient of border forest. Default = 0.11" << endl;
    cout << "-minLeaderLen <value>     Leader length (= height increment) is forced to be at least <value> m" << endl;
    cout << endl;
  }
  //[Usagex]
  ///\endinternal
}//end namespace CrownDensity


int main(int argc, char** argv)
{
  int iterations = 0;//Number of iterations
  double x_coord = 0.0;//Location of the tree
  double y_coord = 0.0;
  string metafile;//File denoting parameters and functions
  int num_parts = 1;    //Number  parts a  segment is  divided  into when
		        //assigning foliage and computing Qabs.
  double hw_start = 0.0;//the age for heartwood formation

  ///\section mainfunction Main growth loop
  //========================================================================
  //Read command line arguments

  if (argc < 3){
    CrownDensity::Usage();
    return 0;
  }
  else{
    iterations = atoi(argv[1]);
    metafile = string(argv[2]);  //< `metafile`, file containing actual parameter files etc.
  }

  //Check possible command line arguments
  //Results will be in an HDF5 file
  if (CheckCommandLine(argc,argv,"-hdf5") == false){
    cout << "No mandatory hdf5 file" <<endl;
    exit(0);
  }
  
  string clarg;
 
  bool toFile = ParseCommandLine(argc,argv,"-toFile", clarg);
  ofstream ff(clarg.c_str() , ofstream::trunc);
  if(!ff)
    toFile = true;   //If opening the output file failed
 
  //The xml file where the tree can be saved
  string xmlfile;
  ParseCommandLine(argc,argv,"-xml",xmlfile);

  //The vertical distribution of fip (from segments)
  string fipfile;
  ParseCommandLine(argc,argv,"-fipdistrib",fipfile);

  //Number parts a segment is  divided into when assigning foliage and
  //computing Qabs.
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-numParts", clarg))
    if(clarg.length() > 0)
      num_parts = atoi(clarg.c_str());
     
  //write interval for printing out result files
  int interval = 0;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-writeInterval", clarg)){
    if (clarg.length()>0){
      interval = atoi(clarg.c_str());
      cout << "Write interval: " << interval <<endl;
    }
  }
 
  //If self_thinning is set, stand density is read from
  //function stems_ha as breast height diameter as argument
  bool self_thinning = false;
  if (CheckCommandLine(argc,argv,"-selfThinning")){
    self_thinning = true;
  }


  //If value of parameter LGAXi is increased during simulation
  string xi;
  bool increase_xi=false;
  int increase_xi_year=0;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-increaseXi",clarg)){
    increase_xi = true;
    increase_xi_year = atoi(clarg.c_str());
  }
  //Implemented  the reading  and use  of tree  files to  generate the
  //forest
  string tree_file;
  if (ParseCommandLine(argc,argv,"-treeFile",tree_file)){
    cout << "Using trees from " <<tree_file <<endl;
  }

  //The heartwood build up starts at age hw_start
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-hw", clarg)) {
    if(clarg.length() > 0) {
      hw_start = atof(clarg.c_str());
      cout  << "HW " << clarg << " "  << hw_start <<endl;
    }
  }

  //Initialize ran3 with some negative integer
  CrownDensity::ran3_seed = -3924678;

  clarg.clear();
  if (ParseCommandLine(argc,argv,"-seed", clarg)){
       CrownDensity::ran3_seed = atoi(clarg.c_str());
      CrownDensity::ran3_seed = -abs(CrownDensity::ran3_seed);
  }
  ran3(&CrownDensity::ran3_seed);


  //Age after which evaluation of tree's self-shading is
  //done with voxelspace calculation (default = never)
  int voxel_calculation = 1000000;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-voxelCalculation",clarg))
    voxel_calculation = atoi(clarg.c_str());

  // If space colonialization is done
  // -space0 = only voxelbox of end of new Segment is checked
  // -space1 = also neighboring boxes in dierction of od the Segments are checked
  // -space2 = also all neighboring boxes are checked
  // Note that space indicators are global variables
  // Note that only one of these should be true at one time
  CrownDensity::space0 = false;
  CrownDensity::space1 = false;
  CrownDensity::space2 = false;
  if(CheckCommandLine(argc,argv,"-space0")) {
    CrownDensity::space0 = true;
  }
  if(CheckCommandLine(argc,argv,"-space1")) {
    CrownDensity::space1 = true;
  }
  if(CheckCommandLine(argc,argv,"-space2")) {
    CrownDensity::space2 = true;
  }
  
  ///\par Parse Extended Borchert-Honda (EBH) allocation of growth
  ///EBH resource distn can be in use in two ways. Both are set by command line
  ///arguments. Option `-EBH` means EBH is in use and values (of lambda parameter)
  ///are specified for all Gravelius orders. Function SPEBHF, in ScotsPine.h,
  ///function file is specified the constructor of the tree. `-EBH1 <value>` means
  ///that EBH is in use and one *value*  is used for all Gravelius orders.
  ///Option `-EBH1 <value>` overrides option `-EBH`. EBH is set by SPis_EBH (Scots Pine Parameter Double SPPD)
  ///to have it in the same way as in CD, 0 == false, 1 == true.
  ///\snippet{lineno} main.cc ebh
  ///\internal
  //[ebh]
  ///\note This condition is checked and set as tree parameter: SetValue(tree, SPis_EBH, value);
  ///\sa Lignum::InitializeTree
  bool growthloop_is_EBH = false;
  if(CheckCommandLine(argc,argv,"-EBH")) {
    growthloop_is_EBH = true;
  }

  bool growthloop_is_EBH1 = false;
  LGMdouble growthloop_EBH1_value = 0.0;
  if (ParseCommandLine(argc,argv,"-EBH1",clarg)) {
    growthloop_is_EBH1 = true;
    growthloop_EBH1_value = atof(clarg.c_str());
    growthloop_is_EBH = true;
  }

  if(growthloop_is_EBH  && growthloop_is_EBH1) { //a bit nonelegant but may work, see initializeTrees()
    ofstream fout("ebh.fun", ofstream::trunc);
    fout << "#Extended Borchert-Honda lambda as a f. of Gravelius order" << endl;
    fout << "#Main axis is 1 first branch 2 etc." << endl;
    LGMdouble x = 1.0;
    for(int i = 0; i < 3; i++) {
      fout << x << " "  <<  growthloop_EBH1_value << endl;
      x += 1.0;
    }
    fout << "6.0 " <<  growthloop_EBH1_value << endl;
    fout.close();
  }
  
  bool growthloop_is_EBH_reduction = false;
  LGMdouble EBH_reduction_parameter = 1.0;
  clarg.clear();
  if(ParseCommandLine(argc,argv,"-EBHREDUCTION", clarg)) {
    growthloop_is_EBH_reduction = true;
    EBH_reduction_parameter = atof(clarg.c_str());
  }

  LGMdouble ebh_final_value = 0.5;
  clarg.clear();
  if(ParseCommandLine(argc,argv,"-EBHFINAL", clarg)) {
    ebh_final_value = atof(clarg.c_str());
  }

  CrownDensity::growthloop_ebh_mode = 0;
  clarg.clear();
  if(ParseCommandLine(argc,argv,"-EBHInput", clarg)) {
    CrownDensity::growthloop_ebh_mode  = atoi(clarg.c_str());
  }
  //[ebh]
  ///\endinternal

  bool growthloop_is_radiation_use_efficiency = false;
  LGMdouble radiation_use_efficiency_parameter = 0.0;
  clarg.clear();
  if(ParseCommandLine(argc,argv,"-RUE", clarg)) {
    growthloop_is_radiation_use_efficiency = true;
    radiation_use_efficiency_parameter = atof(clarg.c_str());
  }


  ParametricCurve af_from_stand;
  ParametricCurve h_from_stand;
  ParametricCurve hcb_from_stand;
  ParametricCurve trees_left_from_stand;
  bool stand_from_file = false;
  if(CheckCommandLine(argc,argv,"-standFromFile")) {
    ifstream input_file("stand-from-file.txt");
    if (!input_file) {
      cout << "Did not find input file stand-from-file.txt" << endl;
      exit(1);
    }

    stand_from_file = true;
    vector<double> age, h, hcb, n, af;
    string line;
    getline(input_file,line);                 //header
    while(input_file.good()){
      getline(input_file,line);
      if(input_file.eof()) {
	break;
      }
      istringstream l(line);
      double ia,ih,ihcb,in,iaf;
      l >> ia >> ih >> ihcb >> in >> iaf;
      age.push_back(ia); h.push_back(ih); hcb.push_back(ihcb);
      n.push_back(in); af.push_back(iaf);
    }  //end of reading input file

    int n_elem = static_cast<int>(age.size());

    af_from_stand = ParametricCurve(n_elem, age, af);
    h_from_stand = ParametricCurve(n_elem, age, h);
    hcb_from_stand = ParametricCurve(n_elem, age, hcb);
    trees_left_from_stand = ParametricCurve(n_elem, age, n);
  }

  CrownDensity::is_adhoc = false;
  if(CheckCommandLine(argc,argv,"-adHoc")) {
    CrownDensity::is_adhoc = true;
  }

  //Use random compoment in segment length
  CrownDensity::is_random_length = false;
  if(CheckCommandLine(argc,argv,"-randomLength")) {
    cout << "Random  component in segment length in use" <<endl;
    CrownDensity::is_random_length = true;
  }

  CrownDensity::is_bud_view_function = false;
  if(CheckCommandLine(argc,argv,"-budViewFunction")) {
    CrownDensity::is_bud_view_function = true;
  }

  clarg.clear();
  if (ParseCommandLine(argc,argv,"-space2Distance",clarg)) {
    CrownDensity::space2_distance = atof(clarg.c_str());
  }


  if(CheckCommandLine(argc,argv,"-heightFun")) {
    CrownDensity::is_height_function = true;
  }

  //bool is_mode_change and int mode_change_yearis are  global variables
  clarg.clear();
  if(ParseCommandLine(argc,argv,"-modeChange", clarg)) {
    CrownDensity::is_mode_change = true;
    CrownDensity::mode_change_year = atoi(clarg.c_str());
  }
  
  int architecture_change_start_year = 0;   //-architectureChange <year> works after this year
  bool architecture_change_is_coming = false;
  clarg.clear();
  //is_architecture_change and architecture_change_year are global variables
  if (ParseCommandLine(argc,argv,"-architectureChange",clarg)){
    Pine::is_architecture_change = false;
    architecture_change_is_coming = true;
    Pine::architecture_change_year = atoi(clarg.c_str());
    architecture_change_start_year = 0;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-aChangeStart",clarg)){
      architecture_change_start_year = atoi(clarg.c_str());
    }   
  }

  bool is_min_leader_length = false;
  double minLeader = 0.0;
  clarg.clear();
  if(ParseCommandLine(argc,argv,"-minLeaderLen", clarg)) {
    is_min_leader_length = true;
    minLeader = atof(clarg.c_str()); 
  }



  
  //See CL argument "-EBH" after instantiatiation of the tree
  //See CL option "-EBH1 <value>" after instantiatiation of the tree

  ///\par  Tree functions 
  ///\snippet{lineno} main.cc TreeInit
  ///\internal
  //[TreeInit]
  //Param sf.fun Specific needle area as function of relative light, domain [0,1] range: meaningfule values for sf  
  //Param fapical.fun Apical dominance as a function of relative light, domain [0,1] in range [0,1]      
  //Param fgo.fun Gravelius order (go) effect on segment length, domain [1, max go] range [0,1]   
  //Param fsapwdown.fun Gravelius order (go) effect passing sapwood down, domain [1, max go] range [0,1]   
  //Param faf.fun parameter estimation? See faf.fun    
  //Param fna.fun Relative light effect on needle angle, domain [0,1] range basically [0,2*pi] but see fna.fun    
  //Param fwd.fun Density of growth ring as function of tree age, domain [0, max age] range [165,214] (based on
  //H. Makinen 2006     
  //Param flr.fun Segment length - radius relationship as function of relative light, domain [0,1] range see flr.fun
  //Param ebh.fun EBH lambda as a function of Gravelius order, domain [1, macx go] range see ebh.fun
  //Param bvf.fun Number of new buds (modifier) as function of local needle area/local volume, domain see bvf.fun
  //range [0,1]
  LignumForest::ScotsPineTree* pine1 = new LignumForest::ScotsPineTree(Point(0.0,0.0,0),PositionVector(0,0,1.0),
								       "sf.fun","fapical.fun","fgo.fun",
								       "fsapwdown.fun","faf.fun","fna.fun", "fwd.fun",
								       "flr.fun", "ebh.fun","bvf.fun");

  //[TreeInit]
  ///\endinternal
  ///\par Bud view function
  ///Bud View Function to Lsystem
  ///Its use is controlled by the global variable is_bud_view_function
  ///If it is == false, bud_view_f is ignored
  ///Bud View Function file is specifed in the constructor of the tree
  ///\snippet{lineno} main.cc BudView
  ///\internal
  //[BudView]
  //Param pine1 The tree
  //Param SPBVF The name of the bud view function. \sa bvf.fun ParametricCurve file 
  CrownDensity::bud_view_f = GetFunction(*pine1, SPBVF);
  //[BudView]
  ///\endinternal
  //SPis_EBH can be set only after the tree has been created
  SetValue(*pine1, SPis_EBH, 0.0);   // i.e. false (1.0 == true)
  if(growthloop_is_EBH) {
    SetValue(*pine1, SPis_EBH, 1.0);
  }

  if (growthloop_is_EBH1) {
    SetValue(*pine1, SPis_EBH, 1.0);
  }

  //Creation of L-system
  Pine::LSystem<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud,PBNAME,PineBudData>* pl1 =
    new Pine::LSystem<LignumForest::ScotsPineSegment,
		      LignumForest::ScotsPineBud,PBNAME,
		      PineBudData>();

  //Heartwood build up age
  SetValue(*pine1,SPHwStart,hw_start);

  // Create an instance of intialization class for tree and initialize the
  // global tree created above
  InitializeTree<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> init_pine1(metafile,VERBOSE);
  init_pine1.initialize(*pine1);

  //The functions for the L-system (they are global variables)

  Pine::fnbuds.install("fnbuds.fun");
  Pine::fnbudslight.install("fnbudslight.fun");
  Pine::fipbud.install("fip-bud.fun");

  
  //Read the amount of incoming radiation (from all directions, no
  // cosine correction) from Firmament instance of the tree and set that
  // as the value of max radiation to the tree
  SetValue(*pine1,TreeQinMax,GetFirmament(*pine1).diffuseBallSensor());

  //Attenuation of radiation inside of shoot (needles) depends
  //on direction expressed by function K(angle)
  ParametricCurve K("K.fun");

  //Function stems_ha specifis the density of the stand as a
  //function of time if self_thinning is not set
  //If self_thinning is set, stems_ha sepcifies density as a
  //function of diameter at breast height
  ParametricCurve stems_ha("stemsha.fun");
 

  //Expand axiom, first initial structure
  pl1->start();
  cout << "Start end: "  << endl;

  //Update Lignum tree pine1 to correspond the structure of the L-system tree
  pl1->lstringToLignum(*pine1,1,PBDATA);

  double wf = 0.0;
  cout << "Collect new foliage begin: " << endl;
  wf = Accumulate(*pine1,wf,CollectNewFoliageMass<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());
  cout << "Collect new foliage end: " << wf << endl;

  //Initial root mass
  SetValue(*pine1,TreeWr,GetValue(*pine1,LGPar)*wf);

  //Calculate  the  LGAsf (specific leaf area) for   newly  created  segments,  sf  in  P
  //Kaitaniemi data depens on segment length
  ForEach(*pine1,SetScotsPineSegmentSf());

  LGMdouble trees_left = stems_ha(0.0);

  //Initialize tree information
  DCLData dcl;
  AccumulateDown(*pine1,dcl,AddBranchWf(),DiameterCrownBase<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());

  cout << "INIT DONE" << endl;

  //======================================================================================
  //   The growth loop
  //=======================================================================================

  double Db_previous = 0.0;
  double Db_current = 0.0;
  ///\par GrowthLoop initialization
  ///GrowthLoop from LignumForest is here to collect data to HDF5
  ///\snippet{lineno} main.cc GLoopInit
  ///\internal
  //[GLoopInit]
  typedef LignumForest::GrowthLoop<LignumForest::ScotsPineTree,LignumForest::ScotsPineSegment,
				   LignumForest::ScotsPineBud,
				   Pine::LSystem<LignumForest::ScotsPineSegment,
						 LignumForest::ScotsPineBud,PBNAME,PineBudData> > ScotsPineGrowthLoop;
  ScotsPineGrowthLoop gloop(pine1,pl1);
  gloop.parseCommandLine(argc,argv);
  gloop.resizeTreeDataMatrix();
  //Initial data
  gloop.collectDataAfterGrowth(0,false);
  //[GLoopInit]
  ///\endinternal
  ///\par HDF5 file initialization
  ///Collect to HDF5 file command line, parameters used, tree functions, all functions, simulation results and  trees in xml format.
  ///Trees must be collected during the simulation to an HDF5 file. The results collected, most notably the 3D matrix for tree data for various attributes,
  ///are saved after simulation.
  ///\snippet{lineno} main.cc HDF5Init
  ///\internal
  //[HDF5Init]
  string hdf5fname;
  ParseCommandLine(argc,argv,"-hdf5", hdf5fname);
  LGMHDF5File hdf5_trees(LignumForest::TREEXML_PREFIX+hdf5fname);
  hdf5_trees.createGroup(LignumForest::TXMLGROUP);
  //[HDF5Init]
  ///\endinternal
  //GROWTH LOOP BEGINS

  cout << "GROWTH LOOP BEGINS" <<endl;
  for (int iter = 0; iter < iterations; iter++){
    cout << "Iter: " << iter << endl;
    //tree age and height to L-system through these global variables
    Pine::L_age = GetValue(*pine1,LGAage);
    Pine::L_H =  GetValue(*pine1,LGAH);

    if(is_mode_change && (iter == mode_change_year)) {     // Change the mode of tree growth
      cout << "Change of morphological mode mode at L_age " << Pine::L_age << endl;

      SetValue(*pine1, SPis_EBH, 0.0);     // --- now FGO or vigor index
      growthloop_is_EBH_reduction = false;
      CrownDensity::is_adhoc = false;                   //No ad hoc lengthening of shoots at crown base

      //This changes parameter values and functions
      InitializeTree<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> init_pine1("MetaFile1.txt",VERBOSE);
      init_pine1.initialize(*pine1);

      Pine::fnbuds = GetFunction(*pine1,LGMNB);  //This for L-System----------     
    }

    //Pine::is_architecture_change influences how morphological development by the L-system works
    if(architecture_change_is_coming && (iter == architecture_change_start_year)) {
      Pine::is_architecture_change = true;
    }

    if(CrownDensity::is_height_function) {
      if(iter == 0) {
	CrownDensity::dDb = 0.003;    // 3 mm --> length growth about 50*0.003 = 0.15
	Db_previous = GetValue(*pine1, LGADbase);
      } else {
	Db_current = GetValue(*pine1, LGADbase);
	CrownDensity::dDb = Db_current - Db_previous;
	Db_previous = Db_current;
      }
    }

    if (increase_xi && iter > increase_xi_year){
      double xii = GetValue(*pine1, LGPxi);
      xii += 0.1/25.0;
      if(xii > 0.85) xii=0.85;
      cout << "Increasing LGPxi from " << GetValue(*pine1,LGPxi) << " to " << xii <<endl;
      SetValue(*pine1,LGPxi, xii);
    }
  
    //The star  mean for each segment 0.14 . 
    LignumForest::SetStarMean<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> setstar(ParametricCurve(0.14));
    ForEach(*pine1,setstar);
 
    double d13cm = GetValue(*pine1,LGADbh)*100.0;

    if (self_thinning)
      trees_left = stems_ha(d13cm);
    else
      trees_left = stems_ha(GetValue(*pine1,LGAage));

    //================================================================
    //Radiation calculations

    //    LGMdouble k_forest = K(PI_VALUE/2.0);
    LGMdouble k_forest = 0.11;
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-kBorderConifer", clarg)){
      k_forest = atof(clarg.c_str());
    }

    LGMdouble Hcb = 0.0;
    LGMdouble treeAf = 0.0;
    treeAf = Accumulate(*pine1,treeAf,CollectFoliageArea<
			LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());

    if(iter > 0) {
      Hcb = dcl.HCrownBase();   //dcl is updated after growth
    }

    double af_to_radiation = treeAf;
    double h_to_radiation = GetValue(*pine1,LGAH);
    double hcb_to_radiation = Hcb;
    global_hcb = Hcb;
    double trees_left_to_radiation = trees_left;
    if(stand_from_file) {
      af_to_radiation = af_from_stand(Pine::L_age);
      h_to_radiation =  h_from_stand(Pine::L_age);
      hcb_to_radiation = hcb_from_stand(Pine::L_age);
      trees_left_to_radiation = trees_left_from_stand(Pine::L_age);
    }

    if(iter < voxel_calculation) {
      cerr << "Year " << iter << " Pairwise radiation calculation" << endl;
	  
      EvaluateRadiationForCfTreeSegmentForest<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> Rad(K,af_to_radiation,k_forest,
													     h_to_radiation,hcb_to_radiation,
													     trees_left_to_radiation,
													     Point(x_coord,y_coord,0));
      ForEach(*pine1,Rad);
    }
    else {
      cerr << "Year " << iter << " Voxel radiation calculation" << endl;
      VoxelSpace vs(Point(0.0,0.0,0.0),Point(1.0,1.0,1.0),
		    0.2,0.2,0.2,5,5,5,GetFirmament(*pine1));

      BoundingBox bb;
      FindCfBoundingBox<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> fb;
      bb = Accumulate(*pine1, bb, fb);

      Point ll = bb.getMin();
      Point ur = bb.getMax();

      vs.resize(ll, ur);
      vs.reset();
      int num_parts = 5;
      bool wood_voxel = true;
      DumpCfTree(vs, *pine1, num_parts, wood_voxel);

      CrownDensity::EvaluateRadiationSelfVoxel<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> Rad(&vs,af_to_radiation,k_forest,h_to_radiation,
													      hcb_to_radiation,trees_left_to_radiation,
													      Point(x_coord,y_coord,0), K);
      ForEach(*pine1,Rad);
    }
    cerr << "CalculateTreeLight end " << endl;

    // If rue*Qin is used in fip, its max value needs to be set

    LGMdouble ini_maxr = 0.0;
    CrownDensity::max_rueqin = Accumulate(*pine1,ini_maxr, LignumForest::FindMaxRueQin());
 
    //===============================================================================
    // Photosynthesis and allocation

    pine1->photosynthesis();

    pine1->respiration();

    //As in LignumForest collect data before growth after respiration
    gloop.collectDataBeforeGrowth(*pine1,0);

    //To print  consistent P and  respirations we must  collect masses
    //now before senescense and new growth
    double wroot = GetValue(*pine1,TreeWr);  //root mass

    
    //Collect foliage 
    LGMdouble wfoliage = 0.0;
    Accumulate(*pine1,wfoliage,CollectFoliageMass<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());

    //Collect sapwood mass
    LGMdouble wsapwood = 0.0;
    Accumulate(*pine1,wsapwood,CollectSapwoodMass<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());

    //TreeAging takes care of senescence in segments (above ground part)
    ForEach(*pine1,TreeAging<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>()); 
    
    //Root mortality
    SetValue(*pine1,TreeWr, 
	     GetValue(*pine1,TreeWr)-GetValue(*pine1,LGPsr)*GetValue(*pine1,TreeWr));

    //Collect  sapwood after  senescence from  all  segments.  Collect
    //again after  new growth excluding new  segments.  The difference
    //of  the two  will tell  how much  sapwood was  needed  in diameter
    //growth
    LGMdouble ws1 = 0.0;
    ws1 = Accumulate(*pine1,ws1,CollectSapwoodMass<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());

    //As in LignumForest update `ws_after_senescence` vector
    UpdateSapwoodAfterSenescence(gloop,ws1,0);
    
    //Uncomment  ForEach  if you  want  to  follow  for each  segment:
    //Radiation   intercepted    and   absorbed;   photosynthesis   and
    //Respiration;  sapwood and foliage  masses; sapwood  an heartwood
    //areas.    cout  <<   "CheckQinQabsPRWsWfAsAh  begin"   <<  endl;
    //ForEach(*pine1,    CheckQinQabsPRWsWfAsAh<ScotsPineSegment,ScotsPineBud>());
    //cout << "CheckQinQabsPRWfAsAh done" << endl;
    
    //This first derive() creates the new segments, whose lengths will be iterated
    pl1->derive();
    pl1->lstringToLignum(*pine1,1,PBDATA);
    //Pass physiological age from mother buds to newly created segments
    double phys_age = 0.0;
    AccumulateDown(*pine1,phys_age,PassPhysiologicalAge<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());
    //Pass the  qin to newly  created segments and to  the terminating
    //buds. Also set LGAip (qin/TreeQinMax) for the terminating buds
    double qin = 0.0;
    PropagateUp(*pine1,qin,LignumForest::ForwardScotsPineQin());

    //Set the needle angle for newly created segments
    ForEach(*pine1,LignumForest::SetScotsPineSegmentNeedleAngle());

    //The variable apical (now attribute of ScotsPineSegment) is set in new segments
    //before length growth by vigor index
    double o0 = 1.0;
    PropagateUp(*pine1,o0,LignumForest::SetScotsPineSegmentApical());
   
    //Calculate the vigour index
    TreePhysiologyVigourIndex(*pine1);

    LGMdouble qin_max = Accumulate(*pine1,qin_max,GetQinMax<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());
    SetValue(*pine1,TreeQinMax,GetFirmament(*pine1).diffuseBallSensor());

    //Length of path from base of tree to each segment
    LGMdouble plength = 0.0;
    PropagateUp(*pine1,plength,PathLength<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());

    //============================================================================
    // Space colonialization
    //=============================================================================
    if (space0 || space1 || space2) {
      BoundingBox bbs;
      FindCfBoundingBox<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> fbs(true); 
      //Segments with needles considered only
      bbs = Accumulate(*pine1, bbs, fbs);

      Point lls = bbs.getMin();
      Point urs = bbs.getMax();

      space_occupancy.resize(lls+Point(-0.5,-0.5,-0.5), urs+Point(0.5,0.5,0.5));
      DumpTreeOccupy<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> dto(3);     //in 3 parts, only foliage parts
      dto.space = &space_occupancy;

      ForEach(*pine1, dto);
    }

    //==============================================================================
    // Extended Borchert-Honda calculation if it is in use
    //==============================================================================

    if(GetValue(*pine1, SPis_EBH) > 0.0) {

      ParametricCurve lambda_fun = GetFunction(*pine1,SPEBHF);
      if(growthloop_is_EBH_reduction) {
	//After age 20 EBH values for all orders change gradually from the nominal value to ebh_final_value 
	if(Pine::L_age > 20) {
	  double v1 = lambda_fun(1.0);
	  double v2 = lambda_fun(2.0);
	  double v3 = lambda_fun(3.0);
	  double v6 = lambda_fun(6.0);
	  LGMdouble p = GetValue(*pine1, LGPapical);
	  v1 *= EBH_reduction_parameter;
	  v2 *= EBH_reduction_parameter;
	  v3 *= EBH_reduction_parameter;
	  v6 *= EBH_reduction_parameter;

	  // if(EBH_reduction_parameter < 1.0) {
	  //    if(v1 < ebh_final_value) v1 = ebh_final_value;
	  //    if(v2 < ebh_final_value) v2 = ebh_final_value;
	  //    if(v3 < ebh_final_value) v3 = ebh_final_value;
	  //    if(v6 < ebh_final_value) v6 = ebh_final_value;
	  //  } else {
	  //    if(v1 > ebh_final_value) v1 = ebh_final_value;
	  //    if(v2 > ebh_final_value) v2 = ebh_final_value;
	  //    if(v3 > ebh_final_value) v3 = ebh_final_value;
	  //    if(v6 > ebh_final_value) v6 = ebh_final_value;
	  //  }
	  if(EBH_reduction_parameter < 1.0) {
	    if(v1 < 0.52) v1 = 0.52;
	    if(v2 < 0.51) v2 = 0.51;
	    if(v3 < 0.503) v3 = 0.503;
	    if(v6 < 0.505) v6 = 0.505;
	  }

	  vector<double> lfo = lambda_fun.getVector();

	  lfo[1] = v1;
	  lfo[3] = v2;
	  lfo[5] = v3;
	  lfo[7] = v6;

	  lambda_fun = ParametricCurve(lfo);
	}
      }
      
      LignumForest::EBH_basipetal_info EBHbI0, EBHbI1;
      //      EBHbI1 = AccumulateDown(*pine1, EBHbI0, EBH_basipetal(lambda_fun) );
      EBHbI1 = AccumulateDown(*pine1, EBHbI0, LignumForest::EBH_basipetal(lambda_fun, CrownDensity::growthloop_ebh_mode) );

      LignumForest::EBH_acropetal_info EBHaI0(1.0, 1.0/lambda_fun(1.0), 1.0);
      PropagateUp(*pine1, EBHaI0, LignumForest::EBH_acropetal(lambda_fun) );

      LignumForest::MaxEBHResource_info m0, m1;
      m0.my_resource = -R_HUGE;

      m1 = AccumulateDown(*pine1, m0, LignumForest::MaxEBHResource() );

      ForEach(*pine1, LignumForest::NormalizeEBHResource(m1.my_resource) );
    }
    
    //=================================================================
    //Allocate    net    photosynthesis    

    //Initialize calculation of thickness growth induced by adding new shoots.
    double alku = 1.0;    //= Gravelius order of main axis
    PropagateUp(*pine1,alku,LignumForest::SetSapwoodDemandAtJunction());



    // REMOVE THESE WHEN YOU REMOVE fgomode,fipmode FROM LGMGrowthAllocator2 !!!!!!!!!!!!!!!!!!
    //ParametricCurve fipmode = GetFunction(*pine1, LGMIP);
    //ParametricCurve fgomode = Lignum::GetFunction(*pine1, Lignum::LGMGO);

    DiameterGrowthData data;
    LGMGrowthAllocator2<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud,LignumForest::SetScotsPineSegmentLength,
			LignumForest::PartialSapwoodAreaDown,LignumForest::ScotsPineDiameterGrowth2,DiameterGrowthData>
      G(*pine1,data,LignumForest::PartialSapwoodAreaDown(GetFunction(*pine1,SPSD)));   


      //Testing the implementation where the sapwood area is passed down
      //as such  between segments that are  in the same  axis.  Only the
      //segments  of higher gravelius  order require  less sapwood  in a
      //branching point.
      try{
	cout << "P=" << G.getP() << " M=" 
	     << G.getM() << " " << " P-M="<< G.getP() - G.getM() << endl;
	cout << "Bisection begin" <<endl;
	Bisection(0.0,100.0,G,0.01,false); //10 grams (C) accuracy 
	cout << "Bisection end " << "L=" << G.getL()<<endl;
      }
      //G will throw an exception if P < M
      catch(TreeGrowthAllocatorException e){
	cout << "P < M " << e.getP() << " " << e.getM() <<endl;
	cout << "Printing last line of output and ending growth loop" <<endl;
	iter = iterations;
      }
      catch(BisectionBracketException e){
	cout << "Could not bracket " << e.getFa() << " " << e.getFb() << " "  << e.getFbl()  <<endl;
	cout << e.getA() << " "  << e.getB() << " " << e.getBl() <<endl;
	cout << "Exit" << endl;
	exit(0);
      }

      if(is_min_leader_length) {
	Axis<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>& stem =  GetAxis(*pine1);
	TreeSegment<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>* last =
	  GetLastTreeSegment(stem);
	double lnew = GetValue(*last, LGAL);
	if(lnew < minLeader) {
	  SetValue(*last, LGAL, minLeader);
	}
      }
  
    //Calculate  the  LGAsf  for   newly  created  segments,  sf  in  P
    //Kaitaniemi data depens on segment length
    ForEach(*pine1,LignumForest::SetScotsPineSegmentSf());
    bool kill = false;
    PropagateUp(*pine1,kill,KillBudsAfterAllocation<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());

    //Now the lengths of the segments are such that Growth = P - M. Do (that is, store permanently,
    //before it was not) the induced diameter growth
    DiameterGrowthData dgdata;
    AccumulateDown(*pine1,dgdata,LignumForest::PartialSapwoodAreaDown(GetFunction(*pine1,SPSD)),
		   LignumForest::ScotsPineDiameterGrowth2(LGMGROWTH));

    double wfnew = 0.0;
    Accumulate(*pine1,wfnew,CollectNewFoliageMass<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());
    SetValue(*pine1,TreeWr, GetValue(*pine1,TreeWr)+GetValue(*pine1,LGPar)*wfnew);
    //To plot how much is allocated to new growth (iWn), diameter growth (iWo)
    //and to roots (iWr) collect datat
    LGMdouble ws2 = 0.0;
    ws2 = Accumulate(*pine1,ws2,CollectOldSegmentSapwoodMass<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());
    LGMdouble ws3 = ws2-ws1;//diameter growth
    //ws4 together with wfnew is new growth
    LGMdouble ws4 = 0.0;
    ws4 = Accumulate(*pine1,ws4,CollectNewSegmentSapwoodMass<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());
    //This is needed for roots
    LGMdouble wr1 = GetValue(*pine1,LGPar)*wfnew;

    //Update L-string, pass the state of the bud to control branching
    pl1->lignumToLstring(*pine1,1,PBDATA);
    pl1->lstringToLignum(*pine1,1,PBDATA);

    //Before derive  pass the  foliage mass of  the mother  segment to
    //terminating buds

    double wftobuds = 0.0;
    PropagateUp(*pine1,wftobuds,ForwardWf<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());

    //============================================================================
    // Here is calculation of estimation of local needle area density for
    // using in estimation of bud fate
    //=============================================================================
    if (is_bud_view_function) {
      BoundingBox b1;
      bool foliage = true;
      FindCfBoundingBox<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> fb1(foliage);
      b1 = Accumulate(*pine1, b1, fb1);

      Point ll = b1.getMin();
      Point ur = b1.getMax();

      VoxelSpace vs1(Point(0.0,0.0,0.0),Point(1.0,1.0,1.0),
		     0.05,0.05,0.05,5,5,5,GetFirmament(*pine1));
      vs1.resize(ll, ur);
      vs1.reset();
      int num_parts = 5;
      DumpCfTree(vs1, *pine1, num_parts);

      LGMdouble cone_height = 0.5, cone_half_angle =  0.7; //= 40 degrees
      int no_points_on_rim = 12;
        
      LignumForest::SetBudViewFunctor sbvf(&vs1, cone_height, cone_half_angle, no_points_on_rim);
      ForEach(*pine1, sbvf);
    }

    pl1->lignumToLstring(*pine1,1,PBDATA);

    //Create new buds in L-system as the function of the mother segment foliage mass
    pl1->derive();
    pl1->lstringToLignum(*pine1,1,PBDATA);
 
    Axis<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>& stem = GetAxis(*pine1);
    //Collect foliage after growth
    LGMdouble wfaftergrowth = 0.0;
    Accumulate(*pine1,wfaftergrowth,CollectFoliageMass<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());
    //Collect sapwood mass after growth
    LGMdouble wsaftergrowth = 0.0;
    Accumulate(*pine1,wsaftergrowth,CollectSapwoodMass<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());
    list<TreeCompartment<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>*>& ls = GetTreeCompartmentList(stem);
    //Mean branch length
    LignumForest::summing bs;
    bs.d2 = bs.d2l = bs.lsum = 0.0;
    bs.n_br = 0;
    bs = Accumulate(*pine1, bs, LignumForest::Branchmeans() );

    //The Qin at the top
    LGMdouble qintop1 = 0.0;
    GetTopQin<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> getTopQin1;
    qintop1 = accumulate(ls.begin(),ls.end(),qintop1, getTopQin1);

    //Diameter and heigth at the crown base.
    AccumulateDown(*pine1,dcl,AddBranchWf(),DiameterCrownBase<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());

    //Crown volume
    CrownVolume<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> cv(0.30);
    double cvol = cv(*pine1);

    //Main stem volume
    MainAxisVolume<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> mav;
    double stemvol = mav(*pine1);

    //Number of segments
    int nsegment = 0;
    nsegment = AccumulateDown(*pine1,nsegment,
			      CountTreeSegments<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>()); 

   

    pl1->prune(*pine1);

    //Set radiation use efficiency in new segments as a function of shadiness
    //experienced by mother segment
    
    if(growthloop_is_radiation_use_efficiency) {   //rue = 1 == no effect by default
      LignumForest::SetRadiationUseEfficiency<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> set_rue(GetFirmament(*pine1).diffuseBallSensor(),
														 radiation_use_efficiency_parameter);
      LGMdouble initial = 0.0;
      PropagateUp(*pine1,initial,set_rue);
    }
    //As in LignumForest Collect data to HDF5 files after growth 
    gloop.collectDataAfterGrowth(iter+1,false);
    CreateTreeXMLDataSet(gloop,hdf5_trees,LignumForest::TXMLGROUP,gloop.getWriteInterval());
  }  // END OF ITERATION END OF ITERATION END OF ITERATION END OF ITERATION

  
  //Clean up.
  cout << "Growth end" << endl;
  cout << "GROWTH DONE " << "NUMBER OF TREES " << gloop.getNumberOfTrees() << endl;
  //Unlike in LignumForest no need for clean up
  ///\par Create HDF5 data
  ///Save collected year by year tree data to an HDF5 file together with data to repeat simulation
  ///\snippet{lineno} main.cc CreateHDF5
  ///\internal
  //[CreateHDF5]
  ParseCommandLine(argc,argv,"-hdf5", hdf5fname);
  TMatrix3D<double>& hdf5_data = gloop.getHDF5TreeData();
  TMatrix2D<double> hdf5_tree_param_data = gloop.getHDF5TreeParameterData();
  CreateLignumHDF5File(hdf5fname,hdf5_data,hdf5_tree_param_data,argc,argv);
  //[CreateHDF5]
  ///\endinternal
  cout << "DATA SAVED TO HDF5 FILES AND SIMULATION DONE" <<endl;
  //Close result file if was open
  if(toFile)
    ff.close();

  pl1->end();  
 
  //Print  Qin and  (x,y,z)  for each  segment  of age  0,1  or 2.  In
  //addition  to Qin,  this gives  you  the idea  of the  size of  the
  //branches
  ForEach(*pine1,PrintSegmentQin<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>());
  LGMdouble wf1 = 0.0;
  cout << "Wf: " << Accumulate(*pine1,wf1,
			       CollectFoliageMass<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>()) << endl;
  LignumForest::summing bs;
  bs = Accumulate(*pine1, bs, LignumForest::Branchmeans());
  if(bs.d2 > 0.0)
    cout << "Averaged2 branch length, m: " << bs.d2l/bs.d2 << endl;
  else
    cout << "Somethin wrong in d2 branch averaging!" << endl;
  

  if(bs.n_br > 0)
    cout << "Average branch length, m: " << bs.lsum/(double)bs.n_br
	 << "   no. branches: " << bs.n_br << endl;
  else
    cout << "Somethin wrong in d2 branch averaging!" << endl;

  int n = 0;
  cout << "Segments: " << Accumulate(*pine1,n, 
				     CountTreeSegments<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>()) << endl;
  LGMdouble qabs = 0.0;
  cout << "Tree       Qabs: " 
       <<  Accumulate(*pine1,qabs,CollectQabs<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud>())<<endl;
 

  //Write the tree to xml file
  if (!xmlfile.empty()){
    XMLDomTreeWriter<LignumForest::ScotsPineSegment,LignumForest::ScotsPineBud> writer;
    writer.writeTreeToXML(*pine1,xmlfile);
  }

  if (!fipfile.empty()){
    cout << "Printing vertical distribution of fip to " << fipfile << endl;
    ofstream fipstream(fipfile.c_str());
    int age = static_cast<int>(GetValue(*pine1,LGAage));
    pair<vector<pair<double,int> >,double> p;
    //number of vertical divisions is tree age
    p.first.resize(age);
    //height intervals is the mean annual growth
    p.second = GetValue(*pine1,LGAH)/GetValue(*pine1,LGAage);
    Accumulate(*pine1,p,LignumForest::CollectVerticalDistributionOfFip());
    double interval = p.second;
    double height = interval;
    //Title
    fipstream  <<left << setfill(' ')
	       << setw(11) << "Height" << setw(11) << "CumulFip" << setw(11) << "Mean fip" <<endl;
    for (unsigned int i = 0; i < p.first.size(); i++){
      if(p.first[i].second > 0){
	fipstream  << setw(11) << height << setw(11) <<  p.first[i].first 
		   << setw(11) << p.first[i].first/p.first[i].second <<endl;
      }
      else{
	fipstream  << setw(11) << height << setw(11) <<  p.first[i].first 
		   << setw(11) << 0 <<endl;
      }
      height = height + interval;
    }
    cout << "Done " <<endl;
  }


  return 0;
#if defined (__APPLE__) || defined(__MACOSX__)
  if(CheckCommandLine(argc,argv,"-viz")) {
    LGMVisualization viz;
    viz.InitVisualization(argc,argv);
    // textures 512x512
    viz.AddCfTree(*pine1, "pine_texture.png", "needle.tga");
    float th = (float)GetValue(*pine1,LGAH);
    cout << th << endl;
    //viz.ResetCameraPosition(th);
    viz.SetMode(SOLID);
    //viz.ResetCameraPosition(GetValue(*pine1,LGAH));
    viz.StartVisualization();
  }
#endif
 

}
