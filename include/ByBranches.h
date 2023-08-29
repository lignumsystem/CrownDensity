#ifndef BYBRANCHES_H
#define BYBRANCHES_H
#include <Lignum.h>


namespace CrownDensity{
template <class TREE>
void makeTreeOfTopWhorls(TREE& pine, TREE* tree);

template <class TREE>
void allocateByBranches(TREE& pine);

template <class TS, class BUD> 
class CollectMainBranches{
 public:
  list<Axis<TS,BUD>*>& operator ()(list<Axis<TS,BUD>*>& b_list, TreeCompartment<TS,BUD>* tc)const {
    if (Axis<TS,BUD>* ax = dynamic_cast<Axis<TS,BUD>*>(tc)){
      list<TreeCompartment<TS,BUD>*>& cmpls =
        GetTreeCompartmentList(*ax);
      if(cmpls.size() > 1) {     //At least one TreeSegment in the Axis
	TreeSegment<TS, BUD>* fs =  GetFirstTreeSegment(*ax);
	if(!(fs == NULL)) {
	  if((GetValue(*fs, LGAomega) == 2.0) && (GetValue(*fs,LGAage) > 2.0)) {
	    b_list.push_back(ax);
	  }
	}
      }
    }
    return b_list;
  }
};

template <class TS, class BUD>
class SetScotsPineSegmentLengthZero {
 public:
  TreeCompartment<TS,BUD>* operator ()(TreeCompartment<TS,BUD>* tc)const {
      if (TS* ts = dynamic_cast<TS*>(tc)){
	if (GetValue(*ts,LGAage) == 0.0){
	  SetValue(*ts, LGAL, 0.0);
	  SetValue(*ts, LGAR, 0.0);
	  SetValue(*ts, LGAWf, 0.0);
	  SetValue(*ts, LGARf, 0.0);
	  SetValue(*ts, LGARh, 0.0);
	}
      }
    return tc;
  }
};

template <class TS, class BUD>
class CountLivingBuds {
 public:
  int& operator ()(int& count, TreeCompartment<TS,BUD>* tc)const {
      if (BUD* bd = dynamic_cast<BUD*>(tc)){
	if (GetValue(*bd,LGAstate) == ALIVE){
	  count++;
	}
      }
    return count;
  }
};

}//End CrownDensity namespace
#endif
