///\file  branchfunctor.cc
///\brief Collect simple descriptive data from branches.
///
///Collect descriptive data from branches with Gravelius order == 2.
#include <SomeFunctors.h>

namespace CrownDensity{ 
///Branch data for Gravelius order 2 (branches from the main trunk)
summing& Branchmeans::operator()(summing& id,
				 TreeCompartment<ScotsPineSegment, ScotsPineBud>*
				 tc)const {
  TreeSegment<ScotsPineSegment,ScotsPineBud>* fs = NULL;
  if (Axis<ScotsPineSegment, ScotsPineBud>* axis =
      dynamic_cast<Axis<ScotsPineSegment, ScotsPineBud>*>(tc)){
    fs = GetFirstTreeSegment(*axis);
    if (fs != NULL){
      if(GetValue(*fs, LGAomega) == 2) {
	std::list<TreeCompartment<ScotsPineSegment, ScotsPineBud>*>&
	  sl = GetTreeCompartmentList(*axis);
	double lb = 0.0;
	list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*>::iterator I
	  = sl.begin();
	///\par Branch length
	///Branch length as a sum of segment lengths
	///\internal
	///\snippet{lineno} branchfunctor.cc LB
	// [LB]
	while(I != sl.end()) {
	  if (ScotsPineSegment* seg = 
	      dynamic_cast<ScotsPineSegment*>(*I)) {
	    lb += GetValue(*seg, LGAL);
	  }
	  I++;
	}
	// [LB]
	///\endinternal
	///\par Diameter data
	///\internal
	///\snippet{lineno} branchfunctor.cc DD
	//[DD]
	//Diameter of the first segment 
	double d = 2.0 * GetValue(*fs, LGAR);
	//Sum of diameter squared from the first segments ("area")
	id.d2 += d * d;
	//Sum of diameter squared times length ("volume")
	id.d2l += d * d * lb;
	//Sum of length of branches
	id.lsum += lb;
	//Sum of branches Gravelius order 2
	id.n_br++;
	//[DD]
	///\endinternal
      }
    }//if fs != NULL
  }//if (Axis....)
  return id;
}
}//end namespace
