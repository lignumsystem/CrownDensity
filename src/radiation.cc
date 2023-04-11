///\file radiation.cc
///\brief Radiation both \f$Q_{in}\f$ and \f$Q_{abs}\f$ for tree segments and buds with border forest.
#include <RadiationCrownDens.h>

#define HIT_THE_FOLIAGE 1
#define NO_HIT 0
#define HIT_THE_WOOD -1



namespace CrownDensity{
///\brief Evaluate self shading using VoxelSpace.
///
///The method is as in LignumForest with `-dumpSelf` command line option.
///Both \f$Q{in}\f$ and \f$Q_{abs}\f$. The borderforest in incoming radiation like in
///stl-lignum EvaluateRadiationForCfTreeSegmentForest
///This evaluates self shading for one ScotsPineSegment, thus
///ForEach does for the whole tree.
///Applicable only for concrete ScotsPineSegment, ScotsPineBud types.
///\todo Move to stl-lignum and implement as a template.
   TreeCompartment<ScotsPineSegment,ScotsPineBud>*
   EvaluateRadiationSelfVoxel::operator()
     (TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
{
  if (ScotsPineSegment* target_segment = dynamic_cast<ScotsPineSegment*>(tc)){
    SetValue(*target_segment, LGAQin, 0.0);
    SetValue(*target_segment, LGAQabs, 0.0);
    //Radiation  conditions are not  evaluated if the segment  has no
    //foliage (in practice  there would be division by  0 in computing
    //absorbed radiation)
    if (GetValue(*target_segment, LGAWf) < R_EPSILON){
	return tc;
    }

    Tree<ScotsPineSegment,ScotsPineBud>& tt = GetTree(*target_segment);
    FirmamentWithMask& firmament = GetFirmament(tt);
    int number_of_sectors = firmament.numberOfRegions();
    vector<double> radiation_direction(3);
    Point middle = GetMidPoint(*target_segment);

    AccumulateOpticalThickness AOT(voxel_space->getXSideLength());
    vector<LGMdouble> s(number_of_sectors,0.0);

    //For forest effect in radiation
    LGMdouble z = middle.getZ();
    LGMdouble dist = sqrt(pow(middle.getX()-stem_loc.getX(),2.0)+
			  pow(middle.getY()-stem_loc.getY(),2.0));
    ///\section RadCond Radiation conditions for a segment
    for (int i = 0; i < number_of_sectors; i++){
      ///\internal
      ///\subsection Qin Incoming radiation
      ///For each sector in Firmament
      ///\snippet{lineno} radiation.cc ForestRad
      ///Incoming radiation with border forest effect: light beam through homogenous border canopy
      // [ForestRad]
      MJ Iop = firmament.diffuseForestRegionRadiationSum
	(i,z,dist,needle_area,forest_k,tree_height,crownbase_height,
	 radiation_direction,density);
      // [ForestRad]
      ///\endinternal
      ///\internal
      ///\snippet{lineno} radiation.cc MaxTrans
      ///Max transmission in voxel
      // [MaxTrans]
      LGMdouble transmission_voxel = 1.0;
      // [MaxTrans]
      ///\endinternal
      ///\internal
      ///Collect light beam path lengths and calculate `optical_thickness`  in voxels in path to sky sector into `vm` vector
      ///\snippet[lineno] radiation.cc VM
      // [VM]
      vector<VoxelMovement> vm;
      PositionVector dir(radiation_direction);
      ParametricCurve K;
      bool dir_star = false;
      voxel_space->getRoute(vm, middle, dir, K, false,dir_star);
      LGMdouble optical_thickness = accumulate(vm.begin(),vm.end(),0.0,AOT);
      // [VM]
      ///\endinternal
      //Vahenna target_segmentin vaikutus pois
      ///\internal
      ///Optical thickness \f$T\f$ for  the segment voxel \post \f$ 0 \leq T \leq 20\f$
      ///\snippet{lineno} radiation.cc OptThick
      // [OptThick]
      LGMdouble k=0.0;
      if(vm[0].n_segs_real > 0.0)
	//Current voxel effect
	k = max(0.0,-0.014+1.056*vm[0].STAR_mean);
      else
	k = 0.0;
      //Optical thickness in the target voxel without targe segment
      optical_thickness -= k * GetValue(*target_segment,LGAAf) * vm[0].l / voxel_space->getBoxVolume();
      // [OptThick]
      ///\endinternal
      if(optical_thickness < 0.0)
	optical_thickness = 0.0;
      ///\internal
      ///Radiation transmission \f$T\f$ in the target voxel \pre  \f$0 \leq T \leq 20\f$
      ///\snippet{lineno} radiation.cc Transm
      // [Transm]
      if(optical_thickness < 20.0){
	//One sector transmission
	transmission_voxel = exp(-optical_thickness);
      }
      else{
	transmission_voxel = 0.0;
      }
      //Transmission, \f$Q_{in}\f$, in one sector
      Iop *= transmission_voxel;
      s[i] = Iop;
      // [Transm]
      ///\endinternal
    } //End of no_sectors ...

    MJ Q_in = accumulate(s.begin(),s.end(),0.0);
    ///\subsection Qabs Absorbed radiation
    ///`s` contains now incoming radiation from each sector. Evaluate how
    ///much segment absorbs from incoming radation.
    LGMdouble Lk, inclination, Rfk, Ack, extinction, sfk, Ask, Wfk;
    Lk = Rfk = Ack =  extinction = sfk = Ask = Wfk = 0.0;
    Lk = GetValue(*target_segment, LGAL);   //length is > 0.0, otherwise we would not bee here
    Rfk = GetValue(*target_segment, LGARf);  //Radius to foliage limit 
    Wfk = GetValue(*target_segment, LGAWf); //Foliage mass
    //sfk  = GetValue(tt, LGPsf); //Foliage m2/kg from tree
    sfk  = GetValue(*target_segment, LGAsf); //Foliage m2/kg from segment!!!
    
    for (int i = 0; i < number_of_sectors; i++){
      firmament.diffuseRegionRadiationSum(i,radiation_direction);
      ///\internal
      ///Calculate the angle,`inclination`, the beam forms the shadow for the segment
      ///\snippet{lineno} radiation.cc SAngle
      // [SAngle]
      LGMdouble a_dot_b = Dot(GetDirection(*target_segment), PositionVector(radiation_direction));
      inclination = PI_DIV_2 - acos(fabs(a_dot_b));
      // [SAngle]
      ///\endinternal
      ///\internal
      ///The shadow area of the segment
      ///\snippet{lineno} radiation.cc SArea
      // [SArea]
      Ack = 2.0*Lk*Rfk*cos(inclination) + PI_VALUE*pow(Rfk,2.0)*sin(inclination);
      // [SArea]
      ///\endinternal
      ///\internal
      ///The empirical extinction coeffient as function of `inclination`
      ///\snippet{lineno} radiation.cc Kincl
      // [Kincl]
      extinction = (double)K(inclination);
      // [Kincl]
      ///\endinternal
      if (Ack == 0.0){
	cout << "ERROR EvaluateRadiationForCfTreeSegment: Ack == 0 (division by 0)"
	     << endl;
      }


      ///\internal
      ///The \f$Q_{abs}\f$  for each sector.
      ///Implement I(k)p = Ip*Ask, Note  Ack must be greater than 0.
      ///It should if there is any foliage.
      ///\snippet{lineno} radiation.cc Ask
      // [Ask]
      Ask = (1.0 - exp(-extinction*((sfk*Wfk)/Ack)))*Ack;
      //QAbs=Qin*Ask
      s[i] *= Ask;
      // [Ask]
      ///\endinternal
    }
    ///\internal
    ///The \f$Q_{abs}\f$ for the segment as sum of all sectors
    ///\snippet{lineno} radiation.cc Qabs
    // [Qabs]
    MJ Q_abs = accumulate(s.begin(),s.end(),0.0);
    SetValue(*target_segment, LGAQabs, Q_abs);
    SetValue(*target_segment, LGAQin, Q_in);
    // [Qabs]
    ///\endinternal
  }
  return tc;
}
}

#undef HIT_THE_FOLIAGE
#undef NO_HIT
#undef HIT_THE_WOOD
