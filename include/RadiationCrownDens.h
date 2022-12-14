#ifndef RADIATIONCROWNDENS_H
#define RADIATIONCROWNDENS_H
#include <Lignum.h>
#include <ScotsPine.h>
#include <VoxelSpace.h>

class EvaluateRadiationSelfVoxel {
public:
  EvaluateRadiationSelfVoxel(VoxelSpace* vs, const LGMdouble& Af, const LGMdouble& for_k,
			     const LGMdouble& tree_h, const LGMdouble& Hcb,
			     const LGMdouble& dens, const Point& x0, const ParametricCurve& k_func):
    voxel_space(vs), needle_area(Af), forest_k(for_k), tree_height(tree_h), crownbase_height(Hcb),
			     density(dens), stem_loc(x0), K(k_func) {}
  TreeCompartment<ScotsPineSegment,ScotsPineBud>*
    operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const;
private:
  VoxelSpace* voxel_space;
  LGMdouble needle_area;
  LGMdouble forest_k;
  LGMdouble tree_height;
  LGMdouble crownbase_height;
  LGMdouble density;
  Point stem_loc;
  ParametricCurve K;
};

class AccumulateOpticalThickness{
 public:
 AccumulateOpticalThickness(LGMdouble side) :
  box_side_length(side) {box_volume = pow(box_side_length,3.0);}
  double operator()(double o_d,VoxelMovement& vm){
    if(vm.af > R_EPSILON && vm.n_segs_real > 0.0) {
      //      LGMdouble k = pow(box_side_length,par_a)*pow(vm.n_segs_real,par_b)*
      //      LGMdouble k = par_a*pow((vm.af/box_volume),par_b)*
      //	max(0.0,-0.014+1.056*vm.STAR_mean);     //NOTE: here transformation STAR_eq
      LGMdouble k = max(0.0,-0.014+1.056*vm.STAR_mean);
      // --> STAR; documented in
      //~/Riston-D/E/LIGNUM/Light/summer-09-test/STAR-vs-STAR_eq.pdf
      o_d += k * vm.af * vm.l / box_volume;
      /* cout << "xn yn zn Af l nseg boxs boxv sST_M k " << vm.x << " " << vm.y << " " << vm.z << " " << vm.af */
      /*      << " " <<  vm.l << " " << vm.n_segs_real << " " << box_side_length << " " << box_volume << " " << */
      /*   vm.STAR_mean << " " << k << endl; */
    }
    return o_d;
  } 
 private:
    LGMdouble box_side_length;
    LGMdouble box_volume;
};



#endif
