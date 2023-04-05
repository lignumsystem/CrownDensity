///\file main.cc
///\brief Main growth loop
///Growth loop steps
///+ Initialize global variables
///+ Read command line
///+ Initialize tree, L-system and global variables from command line
///+ Run the simulation
///+ Save simulation data to HDF5 files.
///\sa CrownDensity
///\page runscript Run crowndensity
///Use the following `run-crowndens.sh` script to run crowndensity.
///Edit `crowndens` command line according to your needs 
///\include run-crowndens.sh
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <Lignum.h>
#include <Bisection.h>
//<From ../LignumForest/
#include <GrowthLoop.h> 
#include <Shading.h>
//From LignumForest project
#include <TreeDataAfterGrowth.h>
//XML file 
#include <XMLTree.h>
//Include the implementation of the tree segment and bud
#include <ScotsPine.h>
#if defined (__APPLE__) || defined(__MACOSX__)
#include <VisualFunctor.h>
//Impelements VisualizeLGMTree
#include <GLSettings.h>
#include <OpenGLUnix.h>
#include <LGMVisualization.h>
#endif
#include <SomeFunctors.h>         //< From ../LignumForest/include
#include <DiameterGrowth.h>       //< From ../LignumForest/include
#include <RadiationCrownDens.h>   //< From ../LignumForest/include
#include <Palubicki_functors.h>   //< From ../LignumForest/include
#include <Space.h>                //< From ../LignumForest/include
//Includes all kinds of stuff, turtle graphics etc.
#include <lengine.h>

//and for pine, see also pine9bp.L in lsys.
namespace Pine{
#include <LSystem.h>

}



using namespace LignumForest;
//#include <ByBranches.h>

///\defgroup g1main Global variables Growth related variables
///@{
int ran3_seed;
bool is_by_branches = false;          ///<this is needed in SomeFunctors.h SetScotsPineSegmentApical

ParametricCurve adhoc("adhoc.fun");   ///<increases growth in lower parts of crown
bool is_adhoc = false;

ParametricCurve toptax("toptax.fun"); ///<function to adjust resource distn among branches

int growthloop_ebh_mode = 0;          ///<This global variable is needed e.g. in bybranches.cc

LGMdouble max_rueqin;
///@}
///\defgroup g2main Global variable space occupance
///@{
///Firmament for space occupancy
Firmament dummy_firm;
///VoxelSpace for space occupancy
VoxelSpace space_occupancy(Point(0.0,0.0,0.0),Point(1.0,1.0,1.0),
			   0.1,0.1,0.1,5,5,5,dummy_firm);
///@}

///\defgroup g3main Global variables L system variables
///These global variables have been declared in pine-em98.L and convey 
//tree age and height to L-system
///@{
extern double L_age, L_H;
///@}
///\defgroup g4main Global variables Height growth
///For conveying height of grown base to SetScotsPineSegmentLength
double global_hcb;      
double dDb;
///@}
///\defgroup g5main Global variables Space colonization
///Space colonialization options for SetScotsPineSegmentLength (in ScotsPine.h)
///and other global variables
///@{
bool space0 = false;
bool space1 = false;
bool space2 = false;
double space2_distance = 0.3;
bool is_forced_height = false;
bool is_height_function = false;
double tax_share = 0.3;        //for toptax
///@}
///\defgroup g6main Global variables Bud view function
///This global variable conveys the Bud View Function to Lsystem
///@{
ParametricCurve bud_view_f;
///If it is in use 
bool is_bud_view_function = false;   
///@}
namespace CrownDensity{
///\section main Main program usage
///\snippet{lineno} main.cc Usagex
///\internal
// [Usagex]
void Usage()
{
  cout << "Usage:  ./crowndens  <iter>  <metafile>" <<endl;
  cout << "-numParts <parts> -hw <hw_start>" <<endl;
  cout << "[-toFile <Filename>] [-xml <filename>]" <<endl;
  cout << "[-fipdistrib <filename>] [-writeInterval interval]" << endl;
  cout << "-[-seed <num>] [-selfThinning] [-increaseXi]" <<endl;
  cout << "[-treeFile <filename>] [-voxelCalculation <year>]" << endl;
  cout << "[-space0] [-space1] [-space2] [-standFromFile] [-adHoc]" << endl;
  cout << "[-budViewFunction] [-EBH] -EBH1 <value>]" << endl;
  cout << "[-space2Distance <Value>] [-byBranches] [-forcedHeight]" << endl;
  cout << "[-heightFun]" << endl;
  cout << endl;
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
  cout << "-heightFun         If length of stem apical shoot is derived from relative crown length (params. LGPe1, LGPe2)." << endl;
}
//[Usagex]
///\endinternal


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
    Usage();
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
  if (CheckCommandLine(argc,argv,"-increaseXi"))
    increase_xi = true;

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
  ran3_seed = -3924678;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-seed", clarg)){
    ran3_seed = atoi(clarg.c_str());
    ran3_seed = -abs(ran3_seed);
  }
  ran3(&ran3_seed);


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

  space0 = false;
  space1 = false;
  space2 = false;
  if(CheckCommandLine(argc,argv,"-space0")) {
    space0 = true;
  }
  if(CheckCommandLine(argc,argv,"-space1")) {
    space1 = true;
  }
  if(CheckCommandLine(argc,argv,"-space2")) {
    space2 = true;
  }

  ///\subsection ebhsection Parse Extended Borchert-Honda (EBH) allocation of growth
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

  growthloop_ebh_mode = 0;
  clarg.clear();
  if(ParseCommandLine(argc,argv,"-EBHInput", clarg)) {
    growthloop_ebh_mode  = atoi(clarg.c_str());
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

  is_adhoc = false;
  if(CheckCommandLine(argc,argv,"-adHoc")) {
    is_adhoc = true;
  }

  is_bud_view_function = false;
  if(CheckCommandLine(argc,argv,"-budViewFunction")) {
    is_bud_view_function = true;
  }

  //  bool is_by_branches = false;
  if(CheckCommandLine(argc,argv,"-byBranches")) {
    is_by_branches = true;
  }

  clarg.clear();
  if (ParseCommandLine(argc,argv,"-space2Distance",clarg)) {
    space2_distance = atof(clarg.c_str());
  }

  if(CheckCommandLine(argc,argv,"-forcedHeight")) {
    is_forced_height = true;
  }

  if(CheckCommandLine(argc,argv,"-heightFun")) {
    is_height_function = true;
  }

  //See CL argument "-EBH" after instantiatiation of the tree

  //See CL option "-EBH1 <value>" after instantiatiation of the tree

  ///\subsection treefunctions  Tree functions 
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
  ScotsPineTree* pine1 = new ScotsPineTree(Point(0.0,0.0,0),PositionVector(0,0,1.0),
					   "sf.fun","fapical.fun","fgo.fun",
					   "fsapwdown.fun","faf.fun","fna.fun", "fwd.fun",
					   "flr.fun", "ebh.fun","bvf.fun");

  //[TreeInit]
  ///\endinternal
  
  ///\subsection budview Bud view function
  ///Bud View Function to Lsystem
  ///Its use is controlled by the global variable is_bud_view_function
  ///If it is == false, bud_view_f is ignored
  ///Bud View Function file is specifed in the constructor of the tree
  ///\snippet{lineno} main.cc BudView
  ///\internal
  //[BudView]
  //Param pine1 The tree
  //Param SPBVF The name of the bud view function. \sa bvf.fun ParametricCurve file 
  bud_view_f = GetFunction(*pine1, SPBVF);
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
  Pine::LSystem<ScotsPineSegment,ScotsPineBud,PineTree::PBNAME,PineTree::PineBudData>* pl1 = new Pine::LSystem<ScotsPineSegment,ScotsPineBud,PineTree::PBNAME,PineTree::PineBudData>();

  //Heartwood build up age
  SetValue(*pine1,SPHwStart,hw_start);

  // Create an instance of intialization class for tree and initialize the
  // global tree created above
  InitializeTree<ScotsPineSegment,ScotsPineBud> init_pine1(metafile,VERBOSE);
  init_pine1.initialize(*pine1);

  
  //    Tämä on hack ---------------------------------------!
  tax_share = GetValue(*pine1, LGPq);   

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
  pl1->lstringToLignum(*pine1,1,PineTree::PBDATA);

  double wf = 0.0;
  cout << "Collect new foliage begin: " << endl;
  wf = Accumulate(*pine1,wf,CollectNewFoliageMass<ScotsPineSegment,ScotsPineBud>());
  cout << "Collect new foliage end: " << wf << endl;

  //Initial root mass
  SetValue(*pine1,TreeWr,GetValue(*pine1,LGPar)*wf);

  //Calculate  the  LGAsf (specific leaf area) for   newly  created  segments,  sf  in  P
  //Kaitaniemi data depens on segment length
  ForEach(*pine1,SetScotsPineSegmentSf());

  LGMdouble trees_left = stems_ha(0.0);

  //Initialize tree information
  DCLData dcl;
  AccumulateDown(*pine1,dcl,AddBranchWf(),DiameterCrownBase<ScotsPineSegment,ScotsPineBud>());

  cout << "Init done" << endl;

  //======================================================================================
  //   The growth loop
  //=======================================================================================

  double Db_previous = 0.0;
  double Db_current = 0.0;
  ///\subsection hdf5 GrowthLoop initialization
  ///GrowthLoop from LignumForest is here to collect data to HDF5
  ///\snippet{lineno} main.cc GLoopInit
  ///\internal
  //[GLoopInit]
  typedef GrowthLoop<ScotsPineTree,ScotsPineSegment,ScotsPineBud,Pine::LSystem<ScotsPineSegment,ScotsPineBud,PineTree::PBNAME,PineTree::PineBudData> > ScotsPineGrowthLoop;
  ScotsPineGrowthLoop gloop(pine1,pl1);
  gloop.parseCommandLine(argc,argv);
  gloop.resizeTreeDataMatrix();
  //Initial data
  gloop.collectDataAfterGrowth(0,false);
  //[GLoopInit]
  ///\endinternal
  ///\subsection hdf5file HDF5 file initialization
  ///Collect to HDF5 file command line, parameters used, tree functions, all functions, simulation results and  trees in xml format
  ///\snippet{lineno} main.cc HDF5Init
  ///\internal
  //[HDF5Init]
  string hdf5fname;
  ParseCommandLine(argc,argv,"-hdf5", hdf5fname);
  LGMHDF5File hdf5_file(hdf5fname);
  LGMHDF5File hdf5_trees(TREEXML_PREFIX+hdf5fname);
  hdf5_file.createGroup(PGROUP);
  hdf5_file.createGroup(TFGROUP);
  hdf5_file.createGroup(AFGROUP);
  hdf5_trees.createGroup(TXMLGROUP);
  ///\sa cxxadt::LGMHDF5File and  TreeDataAfterGrowth.h in LignumForest for HDF5 file group names
  //[HDF5Init]
  ///\endinternal
  for (int iter = 0; iter < iterations; iter++)
  {
    cout << "Iter: " << iter << endl;
    //tree age and height to L-system through these global variables
    L_age = GetValue(*pine1,LGAage);
    L_H =  GetValue(*pine1,LGAH);
       if(iter == 20) {
	 is_adhoc = true;
	 //	 SetValue(pine1, LGPLmin, 0.04);
	 	 is_bud_view_function = true;
		 //	 	 space0 = true;
        }
    // ParametricCurve fgo1;
    // if(iter == 20) {
    //   ParametricCurve apu("fgo1.fun");
    //   fgo1 = apu;
    //   SetFunction(pine1, fgo1, SPFGO);
    // }
    // ParametricCurve fip1;
    // if(iter == 20) {
    //   ParametricCurve apu("fip1.fun");
    //   fip1 = apu;
    //   SetFunction(pine1, fip1, LGMIP);
    // }
    

    //Tree age and height to L-system through these global variables


    if(is_height_function) {
      if(iter == 0) {
	dDb = 0.003;    // 3 mm --> length growth about 50*0.003 = 0.15
	Db_previous = GetValue(*pine1, LGADbase);
      } else {
	Db_current = GetValue(*pine1, LGADbase);
	dDb = Db_current - Db_previous;
	Db_previous = Db_current;
      }
    }

    if (increase_xi){
      if(L_age > 15.0) {
	double xii = GetValue(*pine1, LGPxi);
	xii += 0.1/25.0;
	if(xii > 0.85) xii=0.85;
	cout << "Increasing LGPxi from " << GetValue(*pine1,LGPxi) << " to " << xii <<endl;
	SetValue(*pine1,LGPxi, xii);
      }
    }
  
    //The star  mean for each segment 0.14 . 
    SetStarMean<ScotsPineSegment,ScotsPineBud> setstar(ParametricCurve(0.14));
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

    LGMdouble Hcb = 0.0;
    LGMdouble treeAf = 0.0;
    treeAf = Accumulate(*pine1,treeAf,CollectFoliageArea<
			ScotsPineSegment,ScotsPineBud>());

    if(iter > 0) {
      Hcb = dcl.HCrownBase();   //dcl is updated after growth
    }

      double af_to_radiation = treeAf;
      double h_to_radiation = GetValue(*pine1,LGAH);
      double hcb_to_radiation = Hcb;
      global_hcb = Hcb;
      double trees_left_to_radiation = trees_left;
      if(stand_from_file) {
	af_to_radiation = af_from_stand(L_age);
	h_to_radiation =  h_from_stand(L_age);
	hcb_to_radiation = hcb_from_stand(L_age);
	trees_left_to_radiation = trees_left_from_stand(L_age);
      }

    if(iter < voxel_calculation) {
      cerr << "Year " << iter << " Pairwise radiation calculation" << endl;
	  
      EvaluateRadiationForCfTreeSegmentForest<ScotsPineSegment,ScotsPineBud>
	Rad(K,af_to_radiation,k_forest,h_to_radiation,hcb_to_radiation,
	    trees_left_to_radiation, Point(x_coord,y_coord,0));
      ForEach(*pine1,Rad);
    }
    else {
      cerr << "Year " << iter << " Voxel radiation calculation" << endl;
      VoxelSpace vs(Point(0.0,0.0,0.0),Point(1.0,1.0,1.0),
		    0.2,0.2,0.2,5,5,5,GetFirmament(*pine1));

      BoundingBox bb;
      FindCfBoundingBox<ScotsPineSegment,ScotsPineBud> fb;
      bb = Accumulate(*pine1, bb, fb);

      Point ll = bb.getMin();
      Point ur = bb.getMax();

      vs.resize(ll, ur);
      vs.reset();
      int num_parts = 5;
      bool wood_voxel = true;
      DumpCfTree(vs, *pine1, num_parts, wood_voxel);

      EvaluateRadiationSelfVoxel Rad(&vs,af_to_radiation,k_forest,h_to_radiation,hcb_to_radiation,
				     trees_left_to_radiation, Point(x_coord,y_coord,0), K);
      ForEach(*pine1,Rad);
    }
    cerr << "CalculateTreeLight end " << endl;

    // If rue*Qin is used in fip, its max value needs to be set

    LGMdouble ini_maxr = 0.0;
    max_rueqin = Accumulate(*pine1,ini_maxr, FindMaxRueQin());
 
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
    Accumulate(*pine1,wfoliage,CollectFoliageMass<ScotsPineSegment,ScotsPineBud>());

    //Collect sapwood mass
    LGMdouble wsapwood = 0.0;
    Accumulate(*pine1,wsapwood,CollectSapwoodMass<ScotsPineSegment,ScotsPineBud>());

    //TreeAging takes care of senescence in segments (above ground part)
    ForEach(*pine1,TreeAging<ScotsPineSegment,ScotsPineBud>()); 
    
    //Root mortality
    SetValue(*pine1,TreeWr, 
	     GetValue(*pine1,TreeWr)-GetValue(*pine1,LGPsr)*GetValue(*pine1,TreeWr));

    //Collect  sapwood after  senescence from  all  segments.  Collect
    //again after  new growth excluding new  segments.  The difference
    //of  the two  will tell  how much  sapwood was  needed  in diameter
    //growth
    LGMdouble ws1 = 0.0;
    ws1 = Accumulate(*pine1,ws1,CollectSapwoodMass<ScotsPineSegment,ScotsPineBud>());

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
    pl1->lstringToLignum(*pine1,1,PineTree::PBDATA);

    //Pass the  qin to newly  created segments and to  the terminating
    //buds. Also set LGAip (qin/TreeQinMax) for the terminating buds
    double qin = 0.0;
    PropagateUp(*pine1,qin,ForwardScotsPineQin());

    //Set the needle angle for newly created segments
    ForEach(*pine1,SetScotsPineSegmentNeedleAngle());

    //The variable apical (now attribute of ScotsPineSegment) is set in new segments
    //before length growth by vigor index
    double o0 = 1.0;
    PropagateUp(*pine1,o0,SetScotsPineSegmentApical());
   
    //Calculate the vigour index
    TreePhysiologyVigourIndex(*pine1);

    LGMdouble qin_max = Accumulate(*pine1,qin_max,GetQinMax<ScotsPineSegment,ScotsPineBud>());
    SetValue(*pine1,TreeQinMax,GetFirmament(*pine1).diffuseBallSensor());

    //Length of path from base of tree to each segment
    LGMdouble plength = 0.0;
    PropagateUp(*pine1,plength,PineTree::PathLength<ScotsPineSegment,ScotsPineBud>());

    //============================================================================
    // Space colonialization
    //=============================================================================
    if (space0 || space1 || space2) {
       BoundingBox bbs;
      FindCfBoundingBox<ScotsPineSegment,ScotsPineBud> fbs(true); 
      //Segments with needles considered only
      bbs = Accumulate(*pine1, bbs, fbs);

      Point lls = bbs.getMin();
      Point urs = bbs.getMax();

      space_occupancy.resize(lls+Point(-0.5,-0.5,-0.5), urs+Point(0.5,0.5,0.5));
      DumpTreeOccupy<ScotsPineSegment,ScotsPineBud> dto(3);     //in 3 parts, only foliage parts
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
	if(L_age > 20) {
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
      
      EBH_basipetal_info EBHbI0, EBHbI1;
      //      EBHbI1 = AccumulateDown(*pine1, EBHbI0, EBH_basipetal(lambda_fun) );
       EBHbI1 = AccumulateDown(*pine1, EBHbI0, EBH_basipetal(lambda_fun, growthloop_ebh_mode) );

      EBH_acropetal_info EBHaI0(1.0, 1.0/lambda_fun(1.0), 1.0);
      PropagateUp(*pine1, EBHaI0, EBH_acropetal(lambda_fun) );

      MaxEBHResource_info m0, m1;
      m0.my_resource = -R_HUGE;

      m1 = AccumulateDown(*pine1, m0, MaxEBHResource() );

      ForEach(*pine1, NormalizeEBHResource(m1.my_resource) );
    }
    
    //=================================================================
    //Allocate    net    photosynthesis    

    //Initialize calculation of thickness growth induced by adding new shoots.
    double alku = 1.0;    //= Gravelius order of main axis
    PropagateUp(*pine1,alku,SetSapwoodDemandAtJunction());

    DiameterGrowthData data;
    LGMGrowthAllocator2<ScotsPineSegment,ScotsPineBud,SetScotsPineSegmentLengthBasic,
      PartialSapwoodAreaDown,ScotsPineDiameterGrowth2,DiameterGrowthData>
      G(*pine1,data,PartialSapwoodAreaDown(GetFunction(*pine1,SPSD)));   


    if(is_by_branches && iter > 4) {
      //      allocateByBranches(*pine1);

      if(is_height_function) {
	//Here length of leader at tree level, if height function
	Axis<ScotsPineSegment,ScotsPineBud>& stem =  GetAxis(*pine1);
	TreeSegment<ScotsPineSegment,ScotsPineBud>* last = 
	  GetLastTreeSegment(stem);
	double e1 = GetValue(*pine1,LGPe1);
	double e2 = GetValue(*pine1,LGPe2);
	double Lnew = (e1 + e2*global_hcb/L_H) * dDb;
        if(Lnew < 0.1)           //safeguarding against losing top
            Lnew = 0.1;
	SetValue(*last, LGAL, Lnew);
      }
    } 
    else {   //Normal allocation
      //Testing the implementation where the sapwood area is passed down
      //as such  between segments that are  in the same  axis.  Only the
      //segments  of higher gravelius  order require  less sapwood  in a
      //branching point.
      try{
	cout << "P=" << G.getP() << " M=" 
	     << G.getM() << " " << " P-M="<< G.getP() - G.getM() << endl;
	cout << "Bisection begin" <<endl;
	Bisection(0.0,10.0,G,0.01,/*verbose*/false); //10 grams (C) accuracy 
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

      Axis<ScotsPineSegment,ScotsPineBud>& stem =  GetAxis(*pine1);

      TreeSegment<ScotsPineSegment,ScotsPineBud>* last =
	GetLastTreeSegment(stem);
      double lnew = GetValue(*last, LGAL);
      if(lnew < 0.1)
	SetValue(*last, LGAL, 0.1);

    }

    //Calculate  the  LGAsf  for   newly  created  segments,  sf  in  P
    //Kaitaniemi data depens on segment length
    ForEach(*pine1,SetScotsPineSegmentSf());
    bool kill = false;
    PropagateUp(*pine1,kill,PineTree::KillBudsAfterAllocation<ScotsPineSegment,ScotsPineBud>());

    //Now the lengths of the segments are such that Growth = P - M. Do (that is, store permanently,
    //before it was not) the induced diameter growth
    DiameterGrowthData dgdata;
    AccumulateDown(*pine1,dgdata,PartialSapwoodAreaDown(GetFunction(*pine1,SPSD)),
		   ScotsPineDiameterGrowth2(LGMGROWTH));

    double wfnew = 0.0;
    Accumulate(*pine1,wfnew,CollectNewFoliageMass<ScotsPineSegment,ScotsPineBud>());
    SetValue(*pine1,TreeWr, 
	     GetValue(*pine1,TreeWr)+
	     GetValue(*pine1,LGPar)*wfnew);
    //To plot how much is allocated to new growth (iWn), diameter growth (iWo)
    //and to roots (iWr) collect datat
    LGMdouble ws2 = 0.0;
    ws2 = Accumulate(*pine1,ws2,CollectOldSegmentSapwoodMass<ScotsPineSegment,ScotsPineBud>());
    LGMdouble ws3 = ws2-ws1;//diameter growth
    //ws4 together with wfnew is new growth
    LGMdouble ws4 = 0.0;
    ws4 = Accumulate(*pine1,ws4,CollectNewSegmentSapwoodMass<ScotsPineSegment,ScotsPineBud>());
    //This is needed for roots
    LGMdouble wr1 = GetValue(*pine1,LGPar)*wfnew;

    //Update L-string, pass the state of the bud to control branching
    pl1->lignumToLstring(*pine1,1,PineTree::PBDATA);
    pl1->lstringToLignum(*pine1,1,PineTree::PBDATA);

    //Before derive  pass the  foliage mass of  the mother  segment to
    //terminating buds

    double wftobuds = 0.0;
    PropagateUp(*pine1,wftobuds,PineTree::ForwardWf<ScotsPineSegment,ScotsPineBud>());

    //============================================================================
    // Here is calculation of estimation of local needle area density for
    // using in estimation of bud fate
    //=============================================================================
    if (is_bud_view_function) {
      BoundingBox b1;
      bool foliage = true;
      FindCfBoundingBox<ScotsPineSegment,ScotsPineBud> fb1(foliage);
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
        
      SetBudViewFunctor sbvf(&vs1, cone_height, cone_half_angle, no_points_on_rim);
      ForEach(*pine1, sbvf);
    }

    pl1->lignumToLstring(*pine1,1,PineTree::PBDATA);

    //Create new buds in L-system as the function of the mother segment foliage mass
    pl1->derive();
    pl1->lstringToLignum(*pine1,1,PineTree::PBDATA);
 
    Axis<ScotsPineSegment,ScotsPineBud>& stem = GetAxis(*pine1);

    //Collect foliage after growth
    LGMdouble wfaftergrowth = 0.0;
    Accumulate(*pine1,wfaftergrowth,CollectFoliageMass<ScotsPineSegment,ScotsPineBud>());
    //Collect sapwood mass after growth
    LGMdouble wsaftergrowth = 0.0;
    Accumulate(*pine1,wsaftergrowth,CollectSapwoodMass<ScotsPineSegment,ScotsPineBud>());
    list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*>& ls = GetTreeCompartmentList(stem);
    //Mean branch length
    summing bs;
    bs.d2 = bs.d2l = bs.lsum = 0.0;
    bs.n_br = 0;
    bs = Accumulate(*pine1, bs, Branchmeans() );

    //The Qin at the top
    LGMdouble qintop1 = 0.0;
    PineTree::GetTopQin<ScotsPineSegment,ScotsPineBud> getTopQin1;
    qintop1 = accumulate(ls.begin(),ls.end(),qintop1, getTopQin1);

    //Diameter and heigth at the crown base.
    AccumulateDown(*pine1,dcl,AddBranchWf(),DiameterCrownBase<ScotsPineSegment,ScotsPineBud>());

    //Crown volume
    CrownVolume<ScotsPineSegment,ScotsPineBud> cv(0.30);
    double cvol = cv(*pine1);

    //Main stem volume
    MainAxisVolume<ScotsPineSegment,ScotsPineBud> mav;
    double stemvol = mav(*pine1);

    //Number of segments
    int nsegment = 0;
    nsegment = AccumulateDown(*pine1,nsegment,
			       CountTreeSegments<ScotsPineSegment,ScotsPineBud>()); 

    //There are  some gnuplot  scripts that can  plot data  files from
    //this output.   Please just add  additional output to the  end if
    //needed (please do not  'break' those scripts).  Please note: Due
    //to logic  in the  main loop the  production and  respiration are
    //done before  new growth and  senescence. The masses  plotted are
    //after senescence and new growth.
    if(toFile) {
      double qbs = 0.0;
      double qabs= Accumulate(*pine1,qbs,CollectQabs<ScotsPineSegment,
			    ScotsPineBud>());
//       double tla0 = 0.0;
//       treeAf = Accumulate(*pine1,tla0,CollectFoliageArea<
// 				 ScotsPineSegment,ScotsPineBud>());

      double a0 = 0.0;
      double ASeg0 = Accumulate(*pine1,a0,SurfaceAreaOfNewSegments());

      double wood = 0.0;
      Accumulate(*pine1,wood,CollectWoodMass<ScotsPineSegment,ScotsPineBud>());

      //Collect wood mass in the main axis
      double wstem = GetValue(*pine1,LGAWstem);

      //Wood mass in branches
      double wbranches = wood - wstem;

      //Sapwood mass of stem
      double stem_sw = 0.0;
      Accumulate(*pine1,stem_sw,CollectStemSapwoodMass<ScotsPineSegment,ScotsPineBud>());


      ff <</*1 */ left << setw(6) << setfill(' ') << iter+1 << " "
	 <</*2 */ left << setw(8) << trees_left << " " 
	 <</*3 */ left << setw(8) << GetValue(*pine1,LGAH) << " " 
	 <</*4 */ setw(9) << qintop1 << " " 
	 <</*5 */ setw(12) << qintop1/GetValue(*pine1,TreeQinMax) << " " 
	 <</*6 */ setw(11) << GetValue(*pine1,LGADbase) << " "
	 <</*7 */ setw(11) << GetValue(*pine1,LGADbh) << " "
	 <</*8 */ setw(11) << dcl.DCrownBase() << " " 
	 <</*9 */ setw(11) << dcl.HCrownBase() << " " 
	 <</*10*/ setw(11) << G.getP() << " "      //Production before new growth
	 <</*11*/ setw(11) << G.getM() << " "      //Respiration before new growth: foliage+sapwood+roots
	 <</*12*/ setw(12) << wfaftergrowth << " " //Foliage mass after new growth
	 <</*13*/ setw(14) << bs.d2l/bs.d2 << " "  //Branch means
	 <</*14*/ setw(14) << bs.lsum/(double)bs.n_br  << " "
	 <</*15*/ setw(11) << GetValue(*pine1,TreeWr) << " "     //Root mass after new growth
	 <</*16*/ setw(11) << GetValue(*pine1,LGPmr)*wroot << " "//Root respiration before new growth
	 <</*17*/ setw(11) << G.getM() - GetValue(*pine1,LGPmr)*wroot << " "//Above ground respiration 
	 <</*18*/ setw(11) << wsaftergrowth  << " "           //Sapwood mass after new growth
	 <</*19*/ setw(11) << GetValue(*pine1,LGPms)*wsapwood << " " //Sapwood respiration before new growth
	 <</*20*/ setw(11) << GetValue(*pine1,LGPmf)*wfoliage << " " //Foliage respiration before new growth
	 <</*21*/ setw(11) << cvol << " "  // Crown volume
	 <</*22*/ setw(11) << treeAf << " "//Sum[(LeafA+NeedleA)/Vbox]/Nboxes
	 <</*23*/ setw(11) << nsegment << " " //Number of segments
	 <</*24*/ setw(11) << qabs << " "     // Qabs  
	 <</*25*/ setw(11) << qabs/(GetFirmament(*pine1).
				    diffuseBallSensor()*treeAf) << " "  //Rad eff.
	 <</*26*/ setw(11) << wfoliage << " "  //Foliage that photosynthesized
	 <</*27*/ setw(11) << G.getP()/wfoliage << " "  //The ubiquitous P/Wf
         <</*28*/ setw(11) << treeAf*trees_left << " " //Leaf Area Index
	 <</*29*/ setw(11) << stemvol << " "  //Main axis volume
	 <</*30*/ setw(11) << ASeg0 << " "  //Main axis volume
	 <</*31*/ setw(11) << wfnew << ""   //Foliage produced by elongation (iWn)
	 <</*32*/ setw(11) << ws4 << " "  //Sapwood produced by elongation (iWn)
	 <</*33*/ setw(11) << ws3 << " "  //Sapwood required for diameter growth (iWo)
	 <</*34*/ setw(11) << wr1 << " "  //New roots required by new foliage (iWr)
	 <</*35*/ setw(11) << ws4+ws3 << " " //Sapwood in growth
	 <</*36*/ setw(11) << wood << " " //Wood mass in above-ground parts
	 <</*37*/ setw(11) << wstem << " " //Wood mass of stem
	 <</*38*/ setw(11) << wbranches << " " //Wood mass of branches
	 <</*39*/ setw(11) << stem_sw << " " //Sapwood mass of stem
	 <</*40*/ setw(11) << GetValue(*pine1,LGAAsbase) << " "//Sapwood at base
	 <</*41*/ setw(11) << GetValue(*pine1,LGAAsDbh) << " " //Sapwood at D 1.3
	 <</*42*/ setw(11) << dcl.ASwCrownBase() << " "       //Sapwood at crown base 
	 <</*43*/ setw(11) << G.getL()  //Lambda s.t. G(L) = 0.
	 << endl; 
    }

    bool crowninformation = false;
    if(crowninformation) {
      //Collect and print the two lists of foliage masses and their heights 
      //in the main axis.  Print and assess the crown base afterwards.
      CrownLimitData cld;
      AccumulateDown(*pine1,cld,AddCrownLimitData(),
		     CollectCrownLimitData<ScotsPineSegment,ScotsPineBud>());
      //Collect diameters from the segments in the main axis
      list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*> & tc_ls = GetTreeCompartmentList(GetAxis(*pine1));
      list<double> d_ls;//list of diameters
      d_ls = accumulate(tc_ls.begin(),tc_ls.end(),d_ls,CollectSegmentDiameters());
    
      ostringstream crown_limit_file;
      crown_limit_file << "CrownLimit-" << x_coord << "-" << y_coord << "-" << iter+1 << ".txt";
      //File is CrownLimit+iter+.txt, e.g., "CrownLimit-25.5-22.0-10.txt"
      ofstream cl_file(crown_limit_file.str().c_str());
      const list<pair<double,double> >& hwf_ls = cld.WfHList();
      const list<pair<double,double> >& dh_dwf_ls = cld.dHdWfList();
      const list<pair<double,double> >& h_qabs_ls = cld.HQabsList();
      const list<pair<double,double> >& dh_dqabs_ls = cld.dHdQabsList();
      //How to do this with copy algorithm as with 'cout'??
      list<pair<double,double> >::const_iterator hwf_it = hwf_ls.begin();
      list<double>::const_iterator d_ls_it = d_ls.begin();
      list<pair<double,double> >::const_iterator dh_dwf_it = dh_dwf_ls.begin();
      list<pair<double,double> >::const_iterator h_qabs_it =  h_qabs_ls.begin();
      list<pair<double,double> >::const_iterator dh_dqabs_it = dh_dqabs_ls.begin();

      cl_file << setfill(' ') 
	      << setw(11) << "H" << " " //height
	      << setw(11) << "D" << " " //diameter
	      << setw(11) << "Wf" << " "//foliage
	      << setw(11) << "dH" << " "//growth increment 
	      << setw(11) << "dWf/dH" << " "   
	      << setw(11) << "H" << " " 
	      << setw(11) << "Qabs" << " " 
	      << setw(11) << "dH" << " " 
	      << setw(11) << "dQabs/dH" << " "
	      << setw(11) << "rel H" << endl;
      //For relative height
      LGMdouble Htree = GetValue(*pine1, LGAH);
      while (hwf_it != hwf_ls.end()){
	cl_file << setfill(' ')
		<< setw(11) << (*hwf_it).first << " " 
		<< setw(11) << *d_ls_it << " " 
		<< setw(11) << (*hwf_it).second << " "
		<< setw(11) << (*dh_dwf_it).first << " "
		<< setw(11) << (*dh_dwf_it).second << " "
		<< setw(11) << (*h_qabs_it).first << " "
		<< setw(11) << (*h_qabs_it).second << " "
		<< setw(11) << (*dh_dqabs_it).first << " "
		<< setw(11) << (*dh_dqabs_it).second << " "
		<< setw(11) << (*hwf_it).first / Htree << endl;
	hwf_it++;//the heights and Wf's in the  main axis by branching points 
	dh_dwf_it++;//the same but  now dH and dWf/dH (dWf  is Wf and dH
	//the segment length)
	h_qabs_it++;//the heights and Qabs's the  in the  main axis by branching points 
	dh_dqabs_it++;//the same but now dH and dQabs/dH (dQabs is Qabs and dH
	//the segment length)
	d_ls_it++;
      }

      if (interval && !xmlfile.empty() &&
	  (static_cast<int>(GetValue(*pine1,LGAage)) % interval == 0 )){
	ostringstream xml_interval;
	xml_interval << GetValue(*pine1,LGAage) << "-" << xmlfile;
	XMLDomTreeWriter<ScotsPineSegment,ScotsPineBud> writer;
	cout << "Saving tree to "<< xml_interval.str() << " begin" <<endl; 
	writer.writeTreeToXML(*pine1,xml_interval.str());
	cout << "End" <<endl;
      }
      if (interval && !fipfile.empty()){
	if (static_cast<int>(GetValue(*pine1,LGAage)) % interval == 0){
	  ostringstream fip_interval;
	  fip_interval << GetValue(*pine1,LGAage) << "-" << fipfile << ".txt";
	  cout << "Printing vertical distribution of fip to " << fip_interval.str() << endl;
	  ofstream fipstream(fip_interval.str().c_str());
	  int age = static_cast<int>(GetValue(*pine1,LGAage));
	  pair<vector<pair<double,int> >,double> p;
	  //number of vertical divisions is tree age
	  p.first.resize(age);
	  //height intervals is the mean annual growth
	  p.second = GetValue(*pine1,LGAH)/GetValue(*pine1,LGAage);
	  Accumulate(*pine1,p,CollectVerticalDistributionOfFip());
	  double interval = p.second;
	  double height = interval;
	  //Title
	  fipstream  <<left << setfill(' ')
		     << setw(11) << "Height" << setw(11) << "CumulFip" << setw(11) << "Mean fip" <<endl;
	  for (unsigned int i = 0; i < p.first.size(); i++){
	    if (p.first[i].second > 0){
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
      }
    
    }  //if (crowninformation) {  

    pl1->prune(*pine1);

    //Set radiation use efficiency in new segments as a function of shadiness
    //experienced by mother segment
    
   if(growthloop_is_radiation_use_efficiency) {   //rue = 1 == no effect by default
    SetRadiationUseEfficiency<ScotsPineSegment,ScotsPineBud>
      set_rue(GetFirmament(*pine1).diffuseBallSensor(),
					       radiation_use_efficiency_parameter);
  
    LGMdouble initial = 0.0;
    PropagateUp(*pine1,initial,set_rue);
   }
   //As in LignumForest Collect data to HDF5 files after growth 
   gloop.collectDataAfterGrowth(iter+1,false);
   CreateTreeXMLDataSet(gloop,hdf5_trees,TXMLGROUP,gloop.getWriteInterval());
  }   // END OF ITERATION END OF ITERATION END OF ITERATION END OF ITERATION

  
  //Clean up.
  cout << "Growth end" << endl;
  cout << "GROWTH DONE " << "NUMBER OF TREES " << gloop.getNumberOfTrees() << endl;
  //Unlike in LignumForest no need for clean up
  //gloop.cleanUp();
  ///\subsection collectHDF5 Collect HDF5 data
  ///Save collected year by year tree data to an HDF5 file together with data to repeat simulation
  ///\snippet{lineno} main.cc CollectHDF5
  ///\internal
  //[CollectHDF5]
  TMatrix3D<double>& hdf5_data = gloop.getHDF5TreeData();
  hdf5_file.createDataSet(TREE_DATA_DATASET_NAME,hdf5_data.rows(),hdf5_data.cols(),hdf5_data.zdim(),hdf5_data);
  hdf5_file.createColumnNames(TREE_DATA_DATASET_NAME,TREE_DATA_COLUMN_ATTRIBUTE_NAME,TREE_DATA_COLUMN_NAMES);
  //Parameters used  
  TMatrix2D<double> hdf5_tree_param_data = gloop.getHDF5TreeParameterData();
  hdf5_file.createDataSet(PGROUP+TREE_PARAMETER_DATASET_NAME,hdf5_tree_param_data.rows(),hdf5_tree_param_data.cols(),
			  hdf5_tree_param_data);
  hdf5_file.createColumnNames(PGROUP+TREE_PARAMETER_DATASET_NAME,TREE_PARAMETER_ATTRIBUTE_NAME,TREE_PARAMETER_NAMES);
  //Functions known in a tree
  for (unsigned int i=0; i < FN_V.size();i++){ 
    TMatrix2D<double> hdf5_tree_fn_data = gloop.getHDF5TreeFunctionData(FN_V[i]);
    hdf5_file.createDataSet(TFGROUP+FNA_STR[i],hdf5_tree_fn_data.rows(),hdf5_tree_fn_data.cols(),hdf5_tree_fn_data);
    hdf5_file.createColumnNames(TFGROUP+FNA_STR[i],TREE_FN_ATTRIBUTE_NAME,TREE_FN_COLUMN_NAMES);
  }
  //All functions used
  hdf5_file.createFnDataSetsFromDir("*.fun",AFGROUP,TREE_FN_ATTRIBUTE_NAME,TREE_FN_COLUMN_NAMES);
  //Command line
  vector<string> c_vec;
  std::copy( argv, argv+argc,back_inserter(c_vec));
  ostringstream cline;
  copy(c_vec.begin(),c_vec.end(),ostream_iterator<string>(cline, " "));
  hdf5_file.createDataSet(COMMAND_LINE_DATASET_NAME,cline.str());
  hdf5_file.close();
  //[CollectHDF5]
  ///\endinternal
  cout << "DATA SAVED TO HDF5 FILES AND SIMULATION DONE" <<endl;
  //Close result file if was open
  if(toFile)
    ff.close();

  pl1->end();  
 
  //Print  Qin and  (x,y,z)  for each  segment  of age  0,1  or 2.  In
  //addition  to Qin,  this gives  you  the idea  of the  size of  the
  //branches
  ForEach(*pine1,PineTree::PrintSegmentQin<ScotsPineSegment,ScotsPineBud>());
  LGMdouble wf1 = 0.0;
  cout << "Wf: " << Accumulate(*pine1,wf1,
       		       CollectFoliageMass<ScotsPineSegment,ScotsPineBud>()) << endl;
  summing bs;
  bs = Accumulate(*pine1, bs, Branchmeans());
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
				     CountTreeSegments<ScotsPineSegment,
				     ScotsPineBud>()) << endl;
  LGMdouble qabs = 0.0;
  cout << "Tree       Qabs: " 
       <<  Accumulate(*pine1,qabs,CollectQabs<ScotsPineSegment,ScotsPineBud>())<<endl;
 

  //Write the tree to xml file
  if (!xmlfile.empty()){
    XMLDomTreeWriter<ScotsPineSegment,ScotsPineBud> writer;
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
    Accumulate(*pine1,p,CollectVerticalDistributionOfFip());
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

  //Write information about the grown pine to these files
  ForEach(*pine1, BranchInformation("branchinformation.dat"));
  ForEach(*pine1, SegmentProductionBalance("productionbalance.dat"));

  cout << "Value of Xi at the end " << GetValue(*pine1, LGPxi) << endl;

  // Write information about the branches to file
  ofstream f("branchdata.dat" , ofstream::trunc);
  f << "Height order2 segments segs_w_fol" << endl;
  Sum2ndOrderBranchesSegments sum_segments;
  list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*>& ls 
    = GetTreeCompartmentList(GetAxis(*pine1));
  list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*>::iterator i_c;
  for(i_c = ls.begin(); i_c != ls.end(); i_c++) {
    if(BranchingPoint<ScotsPineSegment,ScotsPineBud>* bp =
       dynamic_cast<BranchingPoint<ScotsPineSegment,ScotsPineBud>*>(*i_c)) {
      list<Axis<ScotsPineSegment,ScotsPineBud>*>& ax_lst = GetAxisList(*bp);
      list<Axis<ScotsPineSegment,ScotsPineBud>*>::iterator i_ax;
      for(i_ax = ax_lst.begin(); i_ax != ax_lst.end(); i_ax++) {  //branches (axes) in whorl
	if(GetTreeCompartmentList(**i_ax).size() > 2)  { //at least one segment
	  Point base = GetPoint(*GetFirstTreeCompartment(**i_ax));
	  vector<int> v(3,0);
	  sum_segments(**i_ax,v);
	  f << base.getZ() << " " << v[0] << " " << v[1] << " " << v[2] << endl;
	}
      }
    }
  }    //for(i_c = ls.begin(); ...

  f.close(); 


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

}//namespace CrownDensity
 
