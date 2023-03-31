//About space colonization etc.

#include <ByBranches.h>
#include <Bisection.h>
#include <DiameterGrowth.h>
#include <XMLTree.h>
#include <Palubicki_functors.h>


extern double tax_share;
namespace CrownDensity{
void allocateByBranches(ScotsPineTree& pine) {

  double H = GetValue(pine, LGAH);
  extern double global_hcb;       //declared and defined in main.cc
  double crown_length = H - global_hcb;
  if(crown_length < 0.1)     //safguarding
     crown_length = 0.1;


  // 1) First collection of main branch axes to a list

  list<Axis<ScotsPineSegment,ScotsPineBud>*> main_branches;

  main_branches = Accumulate(pine, main_branches, CollectMainBranches());

  // 2) Create trees from axes and put them to a list

  ScotsPineTree* tree;
  list<ScotsPineTree*> branch_trees;
  list<Axis<ScotsPineSegment,ScotsPineBud>*>::iterator Ibl;

  for(Ibl = main_branches.begin(); Ibl != main_branches.end(); Ibl++) {
    TreeSegment<ScotsPineSegment,ScotsPineBud>* fs = GetFirstTreeSegment(**Ibl);
    Point start = GetPoint(*fs);
    PositionVector direction = GetDirection(*fs);
    tree = new ScotsPineTree(start,direction, "sf.fun","fapical.fun","fgo.fun",
			     "fsapwdown.fun","faf.fun","fna.fun", "fwd.fun",
			     "flr.fun", "ebh.fun","bvf.fun", **Ibl);
    InitializeTree<ScotsPineSegment,ScotsPineBud> init_pine("MetaFile.txt",QUIET);
    init_pine.initialize(*tree);

    branch_trees.push_back(tree);
  }

  // 3) Tree from two top whorls that have photosythesized, otherwise there are not enough
  // material for growth of the top. The tree will contain whorls with ages 2 and years and
  // a new whorl with age 0

  tree = new ScotsPineTree(Point(0,0,0), PositionVector(0,0,1), "sf.fun","fapical.fun",
			   "fgo.fun", "fsapwdown.fun","faf.fun","fna.fun", "fwd.fun",
			   "flr.fun", "ebh.fun","bvf.fun");
  makeTreeOfTopWhorls(pine, tree);

  branch_trees.push_back(tree);


  // 4) Allocation branch by branch

  // 4.1) Preparations
  extern ParametricCurve toptax;  //defined in main.cc for distribution among branches

  int n_branches = static_cast<int>(branch_trees.size());
  vector<double> surplus(n_branches); 
  vector<double> reduction(n_branches, 0.0);
  vector<double> tax(n_branches, 0.0);   //this is taken from P - M of branches and given back
  
  //tax_share on globaali muuttuja

  vector<int>    n_living_buds(n_branches); 
  vector<bool>   in_iteration(n_branches);
  vector<double> limit_reduction(n_branches, 0.0);
  vector<double> branch_height(n_branches, 0.0);
  double grand_total = 0.0;
  double pool = 0.0, pool0 = 0.0;
  double delta_Ws = 0.0;
  double w_one_shoot = 0.003;  //cost of building one shoot L = 0.251, R = 0.0044
  //Wf = 0.0025, Ws = 0.0023

  list<ScotsPineTree*>::iterator Itl;
  int i = 0;
  for(Itl = branch_trees.begin(); Itl != branch_trees.end(); Itl++) {
    double P = 0.0;
    P = Accumulate(**Itl, P, SumTreePhotosynthesis<ScotsPineSegment,ScotsPineBud>());
    SetValue(**Itl, TreeP, P);
    double M = 0.0;
    M = Accumulate(**Itl, M, SumTreeRespiration<ScotsPineSegment,ScotsPineBud>());
    SetValue(**Itl, TreeM, M);

    grand_total +=  P - M;
 
     if(i < (n_branches - 1) ) {
       Axis<ScotsPineSegment,ScotsPineBud>& axis = GetAxis(**Itl);
       TreeSegment<ScotsPineSegment,ScotsPineBud>* ts1 = GetFirstTreeSegment(axis);
       branch_height[i] = GetPoint(*ts1).getZ();
     } else {   //Tree of four top whorls, its height is the one of the leader shoot
       Axis<ScotsPineSegment,ScotsPineBud>& axis = GetAxis(**Itl);
       TreeSegment<ScotsPineSegment,ScotsPineBud>*ts1 = GetLastTreeSegment(axis);
       branch_height[i] = GetPoint(*ts1).getZ();
     }

    int nb = 0;
    n_living_buds[i] = Accumulate(**Itl,nb,CountLivingBuds());
//    cout << "i nlb limit P-M surplus pool0  ";

    if(n_living_buds[i] == 0) {
      in_iteration[i] = false;
      pool0 += P - M;
      surplus[i] = 0.0;
      limit_reduction[i] = P - M;
      tax[i] = 0.0;
//        cout << i << " " << 0 << " " << 0.0 << " " << P - M << " "
//        << surplus[i] << " " << pool0 << " " << branch_height[i] << endl;

    } else {
      in_iteration[i] = true;
      double limit = static_cast<double>(n_living_buds[i]) * w_one_shoot;
      if((P - M) > limit) {
	pool0 += P - M - limit;
	surplus[i] = limit * (1.0 - tax_share);
	tax[i] = limit * tax_share;
	limit_reduction[i] = P - M - limit;
      } else {
	surplus[i] = (P - M) * (1.0 - tax_share);
	tax[i] = (P - M) * tax_share;
      }

//             cout << i << " " << n_living_buds[i] << " " << limit << " " << P - M << " "
//        << surplus[i] << " " << pool0 << " " << branch_height[i] <<endl;
    }
    i++;
  }
  //  cout << "grand_total pool0  pool% " << grand_total << " " << pool0 << " "
  //       << (int)(100*pool0/(pool0+grand_total)) << endl;

   // 4.1.1 Jaetaan ylimaarainen = grand_total - pool0 = se osuus, jota budit eivat jaksa kayttaa, siis esim.
   // sellaisen oksan tuotos, jossa ei ole yhtaan elavaa silmua.
 
  double additional = 0.0;
  int iteration_branches = 0; 
   for(i = 0; i < n_branches; i++) {
     if(in_iteration[i]) {
       iteration_branches++;
       additional += tax[i];
     }
   }
 
   if(iteration_branches == 0) {    //???????????
     return;
   }
 
  vector<double> addition(n_branches, 0.0); 
  double sum = 0.0;
  for(i = 0; i < n_branches; i++) {
    if(in_iteration[i]) {
      double nb = static_cast<double>(n_living_buds[i]);
      double rh = max(branch_height[i] - global_hcb,0.0) / crown_length;
      addition[i] = toptax(rh) * additional * nb;
      sum += toptax(rh) * nb;
    }
  }

  for(i = 0; i < n_branches; i++) {
    if(in_iteration[i]) {
      addition[i] /= sum;
      //      double rh = max(branch_height[i] - global_hcb,0.0) / crown_length;
       //      cout << "i tot RH addition " << i << " " << additional << " " << rh << " " << addition[i] << endl;
    }
  }

  sum = 0.0;
  for(i = 0; i < n_branches; i++) {
    if(in_iteration[i]) {
      sum += addition[i];
    }
  }
  //  cout << "$$$$$$$$$$$$$   addition0 addition1  " << additional << " " << sum << endl;

  //  cout << "iteration_branches  add   sum  " << iteration_branches << " " << additional << " " << sum << endl;

  // 4.2 Iteration to make utilization in growth == grand_total 

  // 4.2.1 Iteration of branches
  double diff = R_HUGE;
  int iteration = 0;
  vector<double> best_reduction(n_branches,0.0);
  int best_iteration = 0;
  double smallest_diff = R_HUGE;

  while (abs(diff) > grand_total * 0.02) {
    iteration++;
    double f_sum = 0.0, w_sum = 0.0, r_sum = 0.0, p_sum = 0.0;
    int i_branch = 0;
    for(Itl = branch_trees.begin(); Itl != branch_trees.end(); Itl++) {
      if(in_iteration[i_branch]) {

	if(GetValue(pine, SPis_EBH) > 0.0) {

	  ParametricCurve lambda_fun = GetFunction(**Itl,SPEBHF);
	  EBH_basipetal_info EBHbI0, EBHbI1;
	  EBHbI1 = AccumulateDown(**Itl, EBHbI0, EBH_basipetal(lambda_fun) );

	  EBH_acropetal_info EBHaI0(1.0, 1.0/lambda_fun(1.0), 1.0);
	  PropagateUp(**Itl, EBHaI0, EBH_acropetal(lambda_fun) );

	  MaxEBHResource_info m0, m1;
	  m0.my_resource = -R_HUGE;

	  m1 = AccumulateDown(**Itl, m0, MaxEBHResource() );

	  ForEach(**Itl, NormalizeEBHResource(m1.my_resource) );
	}

	DiameterGrowthData data;
	LGMGrowthAllocator2<ScotsPineSegment,ScotsPineBud,SetScotsPineSegmentLength,
	  PartialSapwoodAreaDown,ScotsPineDiameterGrowth2,DiameterGrowthData>
	  G(**Itl,data,PartialSapwoodAreaDown(GetFunction(pine,SPSD)),reduction[i_branch]+
	    limit_reduction[i_branch]-addition[i_branch]); 

	try{
	  p_sum += G.getP() - G.getM();

	  Bisection(0.0,10.0,G,0.01,/*verbose*/false); //10 grams (C) accuracy 
	}
	//If could not bracket
	catch(BisectionBracketException e){
	  ForEach(**Itl, SetScotsPineSegmentLengthZero());
	  //	  cout << "Hep1" << endl;
	}
	DiameterGrowthData data2;
	data2 = AccumulateDown(**Itl,data2,PartialSapwoodAreaDown(GetFunction(pine,SPSD)),
			       ScotsPineDiameterGrowth2(LGMALLOCATE));
	f_sum += GetValue(data2,LGAiWf);
	w_sum += GetValue(data2,LGAiWs);
	r_sum += GetValue(data2,LGAiWf)*GetValue(pine,LGPar);
      }
      i_branch++;
    }
    
    // 4.2.2 Tree level evaluation brings in the need in stem growth

    DiameterGrowthData data1;
    data1 = AccumulateDown(pine,data1,PartialSapwoodAreaDown(GetFunction(pine,SPSD)),
			   ScotsPineDiameterGrowth2(LGMALLOCATE));

    diff = grand_total - GetValue(data1,LGAiWf) - GetValue(data1,LGAiWs) - 
      GetValue(pine,LGPar)* GetValue(data1,LGAiWf);

    if(abs(diff) < smallest_diff) {
      smallest_diff = abs(diff);
      for(i = 0; i < n_branches; i++) {
	best_reduction[i] = reduction[i_branch]+limit_reduction[i_branch]-addition[i_branch];
      }
    }
    
//     cout << "**All** P - M  iWs  iWf iWr sum  iteration diff " << grand_total << " "
// 	 << GetValue(data1,LGAiWs) << " "
// 	 <<  GetValue(data1,LGAiWf) <<  " " << GetValue(pine,LGPar)* GetValue(data1,LGAiWf) << " "
// 	 << GetValue(data1,LGAiWf) + GetValue(data1,LGAiWs) + GetValue(pine,LGPar)* GetValue(data1,LGAiWf)
// 	 << " " << iteration << " " << diff << endl;

    double bud_sum = 0.0;
    vector<double> d_reduction(n_branches, 0.0);  
    for(i = 0; i < n_branches; i++) {
      if(in_iteration[i]) {
	double nb = static_cast<double>(n_living_buds[i]);
	if(diff < 0.0) {
	  d_reduction[i] = nb * (-diff/4.0);
	  bud_sum += nb;
	} else {
	  d_reduction[i] = nb * (-diff/4.0);
	  bud_sum += nb;
	}
      }
    }
    for(i = 0; i < n_branches; i++) {
      if(in_iteration[i]) {
	reduction[i] += d_reduction[i]/bud_sum;
      }
    }

    if(iteration > 20)
      break;
  }   //while(diff ...

  if(iteration > 20) {   //max iterations exceeded, go with smallest difference so far
    for(Itl = branch_trees.begin(); Itl != branch_trees.end(); Itl++) {
      if(in_iteration[i]) {
	DiameterGrowthData data;
	LGMGrowthAllocator2<ScotsPineSegment,ScotsPineBud,SetScotsPineSegmentLength,
	  PartialSapwoodAreaDown,ScotsPineDiameterGrowth2,DiameterGrowthData>
	  G(**Itl,data,PartialSapwoodAreaDown(GetFunction(pine,SPSD)),best_reduction[i]); 

	try{
	  Bisection(0.0,10.0,G,0.01,/*verbose*/false); //10 grams (C) accuracy 
	}
	catch(BisectionBracketException e){
	  ForEach(**Itl, SetScotsPineSegmentLengthZero());
	}
      }
    }

    DiameterGrowthData data1;
    data1 = AccumulateDown(pine,data1,PartialSapwoodAreaDown(GetFunction(pine,SPSD)),
			   ScotsPineDiameterGrowth2(LGMALLOCATE));
  }



  return;
}


void makeTreeOfTopWhorls(ScotsPineTree& pine, ScotsPineTree* tree) {

  Axis<ScotsPineSegment,ScotsPineBud>& stem =  GetAxis(pine);
  list<TreeCompartment<ScotsPineSegment, ScotsPineBud>*>& stem_list = 
    GetTreeCompartmentList(stem);
  list<TreeCompartment<ScotsPineSegment, ScotsPineBud>*>::iterator Istem = 
    stem_list.end();

  Axis<ScotsPineSegment,ScotsPineBud>& new_stem =  GetAxis(*tree); 
//   list<TreeCompartment<ScotsPineSegment, ScotsPineBud>*>& new_stem_list = 
//     GetTreeCompartmentList(new_stem);

  // Eight items, that is, Bud, TS, BP, TS, BP, TS, BP, TS are copied from stem_list
  // to new_stem list. NOTE that as the tree (pine) is between the two passes of
  // L-system, there is no BP between Bud and TS; after the second pass there will be
  // a BP between Bud and TS.

  for(int i = 0; i < 8; i++) {
    Istem--;
  }           //*Istem now points to the lowest of 8 TreeCompartments in the list
  for(int i = 0; i < 8; i++) {
    InsertTreeCompartment(new_stem, *Istem);
    Istem++;
  }


  return;
}
}
