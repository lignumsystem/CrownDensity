#ifndef BYBRANCHES_H
#define BYBRANCHES_H
#include <Lignum.h>
#include <ScotsPine.h>

void makeTreeOfTopWhorls(ScotsPineTree& pine, ScotsPineTree* tree);

void allocateByBranches(ScotsPineTree& pine);


class CollectMainBranches {
 public:
  list<Axis<ScotsPineSegment,ScotsPineBud>*>&
    operator ()(list<Axis<ScotsPineSegment,ScotsPineBud>*>& b_list,
		TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const {
    if (Axis<ScotsPineSegment,ScotsPineBud>* ax =
        dynamic_cast<Axis<ScotsPineSegment,ScotsPineBud>*>(tc)){
      list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*>& cmpls =
        GetTreeCompartmentList(*ax);
      if(cmpls.size() > 1) {     //At least one TreeSegment in the Axis
	TreeSegment<ScotsPineSegment, ScotsPineBud>* fs =  GetFirstTreeSegment(*ax);
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


class SetScotsPineSegmentLengthZero {
 public:
  TreeCompartment<ScotsPineSegment,ScotsPineBud>*
    operator ()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const {
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
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

class CountLivingBuds {
 public:
  int& operator ()(int& count, TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const {
      if (ScotsPineBud* bd = dynamic_cast<ScotsPineBud*>(tc)){
	if (GetValue(*bd,LGAstate) == ALIVE){
	  count++;
	}
      }
    return count;
  }
};


#endif
