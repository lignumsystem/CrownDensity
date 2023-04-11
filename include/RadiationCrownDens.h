///\file RadiationCrownDens.h
#ifndef RADIATIONCROWNDENS_H
#define RADIATIONCROWNDENS_H
#include <Lignum.h>
#include <ScotsPine.h>
#include <VoxelSpace.h>

namespace CrownDensity{
///\brief Evaluate self shading using VoxelSpace.
class EvaluateRadiationSelfVoxel {
public:
  ///\param vs Voxel space
  ///\param Af Segment foliage area
  ///\param for_k Forest extinction parameter for K function
  ///\param tree_h Tree height
  ///\param Hcb Height crown base
  ///\param dens Homogenous border forest density
  ///\param x0 Location of the Tree stem
  ///\param k_func Empirical extinction function K (also of `for_k`)
  EvaluateRadiationSelfVoxel(VoxelSpace* vs, const LGMdouble& Af, const LGMdouble& for_k,
			     const LGMdouble& tree_h, const LGMdouble& Hcb,
			     const LGMdouble& dens, const Point& x0, const ParametricCurve& k_func):
    voxel_space(vs), needle_area(Af), forest_k(for_k), tree_height(tree_h), crownbase_height(Hcb),
			     density(dens), stem_loc(x0), K(k_func) {}
  TreeCompartment<ScotsPineSegment,ScotsPineBud>*
    operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const;
private:
  VoxelSpace* voxel_space;///< Voxel space
  LGMdouble needle_area; ///< Segment needle area
  LGMdouble forest_k; ///< Forest extinction coeffient \sa K
  LGMdouble tree_height;///<Tree height
  LGMdouble crownbase_height;///Crown base height
  LGMdouble density; ///Border forest density
  Point stem_loc;///Tree main stem location
  ParametricCurve K;///<Empirical extinction function (also of `forest_k`)\sa forest_k
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

}

#endif
