///\file ScotsPine.cc
///\brief Photosynthesis, respiration and aging
#include <ScotsPine.h>

namespace CrownDensity{
///\brief Photosynthetic production.
///Photosynthetic production implemented as 
///\f[
///P = rue \times LGP_{pr} \times Q_{abs}
///\f]
///where parameter `Lignum::LGPpr` = Photosynthetic efficiency
///(see LGMSymbols.h). `rue` depends on the radiation conditions of the CfTreeSegment
///at its birth. At full light (at top of the stand) `rue` = 1, in shaded
///conditions possibly `rue` > 1.
///\pre `CrownDensity::SPrue` must be defined
///\todo ScotsPineSegment::photosynthesis(), ScotsPineSegment::respiration() and ScotsPineSegment::aging() could be implemented as functors
/// to ease the testing of experimental models.
void ScotsPineSegment::photosynthesis()
{
  Tree<ScotsPineSegment,ScotsPineBud>& t = GetTree(*this);
  SetValue(*this,LGAP, GetValue(*this, SPrue) * GetValue(t, LGPpr) * GetValue(*this,LGAQabs));

  // //const ParametricCurve& fip = GetFunction(GetTree(*this),LGMIP);
  // //  double f_ip = fip(GetValue(*this,LGAQin)/GetValue(t,TreeQinMax));
  // double ip = GetValue(*this,LGAQin)/GetValue(t,TreeQinMax);
  // double pr = GetValue(t, LGPpr);
  // double delta = pr - 0.0006;
  // //delta *= f_ip;
  // double f_ip = 0.0;
  // if(ip < 0.6)
  //   f_ip = 0.0;
  // else
  //   f_ip = (ip - 0.6)/0.4;
  // delta *= f_ip;
  
  // pr = 0.0006 +  delta;

  // SetValue(*this, LGAP, pr*GetValue(*this,LGAQabs));
}

///\brief Respiration
///
///Respiration rate of the segment as the function of needle mass
///and sapwood mass. Respiration is
///\f[
/// R_{tot} = R_{W_f} + R_{W_s}
///\f]
///where
///\f[
///\begin{cases}
/// R_{W_f} = LGP_{mf} \times W_f\\
/// R_{W_s} = LGP_{ms} \times {W_s}.
///\end{cases}
///\f]
///\todo Redefine respiration, that is defined in CfTreeSeegmentMetabolismI.h
///to take into account wood density that varies between annual rings
///according to function `SPWD`.
///The variation of wood density among tree rings is taken into
///consideration in an average manner: mean density of sapwood is
///calculated with the aid of function `SPWD` that defines density of
///annual ring as a function of age of segment (at creation of annual
///ring). The average is taken of function `SPWD` over range (hardwood
//radius/radius) x segment age ... segment age.
void ScotsPineSegment::respiration()
{
  LGMdouble resp = 0.0;
  Tree<ScotsPineSegment,ScotsPineBud>& t = GetTree(*this);
  //Rtot = Rfoliage + Rsapwood
  resp = GetValue(t, LGPmf)*GetValue(*this,LGAWf); //Foliage
  
//   const ParametricCurve& wd = GetFunction(dynamic_cast<ScotsPineTree&>(t), SPWD);
//   LGMdouble age = GetValue(*this, LGAage);
//   LGMdouble rh = GetValue(*this, LGARh);
//   LGMdouble r = GetValue(*this, LGAR);
//   int start = (int)((rh/r)*age);
//   LGMdouble dens = 0.0;
//   LGMdouble count = 0.0;
// 	//This is a shorthand fix: Hakkila 1969 CIFF 67.6 shows that wood density for branches
// 	// is > 200 kgC/m3. To achieve this add 10 years to age, if branch

//   for(int a = start; a <= (int)age; a++) {
//     LGMdouble aa = a;
//     if(GetValue(*this, LGAomega) > 1.0)
//       aa += 10.0;
//     dens += wd((double)aa);
//     count += 1.0;
//   }
//   if(count < 0.9999999)
//     dens = wd(0.0);
//   else
//     dens /= count;

  //resp += GetValue(t,LGPms)*dens*GetValue(*this,LGAAs)*GetValue(*this,LGAL);
  resp += GetValue(t,LGPms)*GetValue(*this,LGAWs);//Sapwood respiration
  SetValue(*this,LGAM, resp);
}

///Aging  tailored for ScotsPineSegment:  heartwood formation  after `CrownDensity::SPHWStart` years.
///\pre Sapwood senescence age `CrownDensity::SPHwStart` must be defined (from command line )
///\pre Foliage senescence function `Lignum::LGMFM` must be defined
///\sa run-crowndens.sh
void ScotsPineSegment::aging()
{
  ///\internal
  ///**The aging steps**
  ///+ Add age 1 year (time step)
  ///+ Sapwood senescence if `LGAage` > `SPHwStart` (e.g. 13 years, Bjorklund, Silva Fennica)
  ///+ Foliage senescence as a declining function of age
  /// \sa CrownDensity::SPHwStart Lignum::LGMFM
  ///\endinternal
  //One year older
  SetValue(*this,LGAage,GetValue(*this,LGAage)+1.0);
  //Sapwood senescence
  if (GetValue(*this,LGAage) > GetValue(dynamic_cast<const ScotsPineTree&>(GetTree(*this)),SPHwStart)){
    LGMdouble dAs = GetValue(GetTree(*this),LGPss) * GetValue(*this,LGAAs);
    LGMdouble Ah_new =  dAs + GetValue(*this, LGAAh);
    LGMdouble Rh_new = sqrt(Ah_new/PI_VALUE);
    SetValue(*this,LGARh,Rh_new);
  }
  //Foliage mass decline
  const ParametricCurve& fm = GetFunction(GetTree(*this),LGMFM);
  //This implementation assumes declining function of age from 1 to 0.
  LGMdouble Wf_new = fm(GetValue(*this,LGAage))*GetValue(*this,LGAWf0);
  SetValue(*this,LGAWf,Wf_new);
  
}
}
