#ifndef CROWNDENSITYGLOBALS_H
#define CROWNDENSITYGLOBALS_H
///\file CrownDensityGlobals.h
///\brief Global varables in CrownDensity
///\sa crowndensity_globals.cc

#include <LGMSymbols.h>
#include <Point.h>
#include <ParametricCurve.h>
#include <Firmament.h>
#include <VoxelSpace.h>

///\defgroup lsysglobal L-system global variables
///@{
///These global variables have been declared in pine-em98.L and convey 
///tree age and height to L-system
extern double L_age;
extern double L_H;
///@}

extern int ran3_seed;

extern bool is_by_branches;          

extern ParametricCurve adhoc;   
extern bool is_adhoc;

extern int growthloop_ebh_mode;       

extern LGMdouble max_rueqin;

extern bool is_mode_change;
extern int mode_change_year;    

extern Firmament dummy_firm;
extern VoxelSpace space_occupancy;

//extern double L_age, L_H;

extern double global_hcb;      
extern double dDb;

extern bool space0;
extern bool space1;
extern bool space2;
extern double space2_distance;
extern bool is_forced_height;
extern bool is_height_function;

extern ParametricCurve toptax; ///<function to adjust resource distn among branches
extern double tax_share;        

extern ParametricCurve bud_view_f;
extern bool is_bud_view_function;   
#endif
