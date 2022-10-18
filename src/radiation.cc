#include <RadiationCrownDens.h>

#define HIT_THE_FOLIAGE 1
#define NO_HIT 0
#define HIT_THE_WOOD -1




//This functor EvaluateRadiationSelfVoxel evaluates self shading
//using VoxelSpace: as in LignumForest with -dumpSelf command line option.
//Both Qin ja Q abs.
//This evaluates self shading for one ScotsPineSegment, thus
//ForEach does for the whole tree
//Now for ScotsPineSegment, ScotsPineBud move maybe later 
//to stl-lignum as a template

//Now the borderforest in incoming radiation like in
//stl-lignum EvaluateRadiationForCfTreeSegmentForest


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


    for (int i = 0; i < number_of_sectors; i++){
      MJ Iop = firmament.diffuseForestRegionRadiationSum
	(i,z,dist,needle_area,forest_k,tree_height,crownbase_height,
	 radiation_direction,density);
  
      LGMdouble transmission_voxel = 1.0;
      vector<VoxelMovement> vm;
      PositionVector dir(radiation_direction);
      ParametricCurve K;
      bool dir_star = false;
      voxel_space->getRoute(vm, middle, dir, K, false,dir_star);
      LGMdouble optical_thickness = accumulate(vm.begin(),vm.end(),0.0,AOT);

      //Vahenna target_segmentin vaikutus pois 
      LGMdouble k;
      if(vm[0].n_segs_real > 0.0)
	k = max(0.0,-0.014+1.056*vm[0].STAR_mean);
      else
	k = 0.0;
      optical_thickness -= k * GetValue(*target_segment,LGAAf) *
	vm[0].l / voxel_space->getBoxVolume();

      if(optical_thickness < 0.0)
	optical_thickness = 0.0;

      if(optical_thickness < 20.0){
	transmission_voxel = exp(-optical_thickness);
      }
      else{
	transmission_voxel = 0.0;
      }

      Iop *= transmission_voxel;
      s[i] = Iop;
    } //End of no_sectors ...

    MJ Q_in = accumulate(s.begin(),s.end(),0.0);
    
    //s contains now incoming radiation from each sector. Evaluate how
    //much segment absorbs from incoming radation.
    LGMdouble Lk, inclination, Rfk, Ack, extinction, sfk, Ask, Wfk;
    Lk = Rfk = Ack =  extinction = sfk = Ask = Wfk = 0.0;
    Lk = GetValue(*target_segment, LGAL);   //length is > 0.0, otherwise we would not bee here
    Rfk = GetValue(*target_segment, LGARf);  //Radius to foliage limit 
    Wfk = GetValue(*target_segment, LGAWf); //Foliage mass
    //sfk  = GetValue(tt, LGPsf); //Foliage m2/kg from tree
    sfk  = GetValue(*target_segment, LGAsf); //Foliage m2/kg from segment!!!

    for (int i = 0; i < number_of_sectors; i++){
      firmament.diffuseRegionRadiationSum(i,radiation_direction);
      LGMdouble a_dot_b = Dot(GetDirection(*target_segment), PositionVector(radiation_direction));
      inclination = PI_DIV_2 - acos(fabs(a_dot_b));

      Ack = 2.0*Lk*Rfk*cos(inclination) + PI_VALUE*pow(Rfk,2.0)*sin(inclination);
      extinction = (double)K(inclination);

      if (Ack == 0.0){
	cout << "ERROR EvaluateRadiationForCfTreeSegment: Ack == 0 (division by 0)"
	     << endl;
      }

      //implement I(k)p = Ip*Ask, Note  Ack must be greater than 0 (it
      //should if there is any foliage)
      Ask = (1.0 - exp(-extinction*((sfk*Wfk)/Ack)))*Ack;
      s[i] *= Ask;
    }
    MJ Q_abs = accumulate(s.begin(),s.end(),0.0);
    SetValue(*target_segment, LGAQabs, Q_abs);
    SetValue(*target_segment, LGAQin, Q_in);
  }
  return tc;
}


#undef HIT_THE_FOLIAGE
#undef NO_HIT
#undef HIT_THE_WOOD
