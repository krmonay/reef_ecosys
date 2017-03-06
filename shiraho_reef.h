/*
** svn $Id: inlet_test.h 585 2012-01-03 18:44:28Z arango $
*******************************************************************************
** Copyright (c) 2002-2012 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Inlet Test Case, waves-ocean (SWAN/ROMS) two-way coupling.
**
** Application flag:   SHIRAHO_REEF
** Input script:       ocean_inlet_test.in
**                     coupling_inlet_test.in
**                     sediment_inlet_test.in
*/

#define UV_VIS2
/*#define MIX_S_UV*/
#define MIX_GEO_UV
#define MIX_GEO_TS
#undef  ANA_GRID
#define MASKING
#define UV_ADV
#define UV_COR

/*#define UV_U3ADV_SPLIT*/
/*#define UV_C2ADVECTION*/
/*#define UV_C4ADVECTION*/
/*#define UV_SADVECTION*/

/*#define TS_U3ADV_SPLIT*/
/*#define TS_A4HADVECTION*/
/*#define TS_C4HADVECTION*/
/*#define TS_MPDATA*/
#define TS_U3HADVECTION
/*#define TS_A4VADVECTION*/
/*#define TS_C4VADVECTION*/
/*#define TS_C2VADVECTION*/
/*#define TS_SADVECTION*/
/*#define TS_SMAGORINSKY*/

#define NONLIN_EOS
#define SALINITY
#define DJ_GRADPS
#define FSOBC_REDUCED
#define SOLVE3D
/*#define SPLINES*/
#define WET_DRY

/*#define AVERAGES*/

#define ANA_INITIAL
/*#define ANA_FSOBC*/
#define ANA_M2OBC
#define ANA_TOBC

#define SOLAR_SOURCE

#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE
# define EMINUSP
/*# define ANA_CLOUD*/
/*# define ANA_HUMID*/
/*# define ANA_PAIR*/
/*# define ANA_TAIR*/
/*# define ANA_RAIN*/
/*# define ANA_WINDS*/
#else
# define ANA_SMFLUX
# define ANA_STFLUX
#endif

#define DIAGNOSTICS_UV
#define DIAGNOSTICS_TS

/* river discharge */

/*#define UV_PSOURCE*/
#ifdef UV_PSOURCE
# define TS_PSOURCE
/*# define ANA_PSOURCE*/
#endif

/* submarine groundwater discharge */

/*#define SGD_ON*/    /*Original CPP flag */

/* waves-ocean (SWAN/ROMS) two-way coupling. */

#define SWAN_COUPLING
#define NEARSHORE_MELLOR08
#ifdef SWAN_COUPLING
# define MCT_LIB
#endif
/*#define AIR_OCEAN*/

/* define only one of the following 5 */
#define UV_LOGDRAG
#undef  UV_QDRAG
#undef  MB_BBL
#undef  SG_BBL
#undef  SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
#endif

#ifdef SOLVE3D
# define GLS_MIXING
/*# define MY25_MIXING*/
# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  undef  CANUTO_A
/*#  define K_C4ADVECTION*/
#  define N2S2_HORAVG
# endif
# define SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  undef  BEDLOAD_SOULSBY
#  undef  BEDLOAD_MPM
#  define SED_MORPH
# endif
# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  define ANA_SEDIMENT
#  define REVER_SEDIMENT
# endif
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
#endif

/***  Biological model options. ***/

#define REEF_ECOSYS  /*Original CPP flag */
#ifdef REEF_ECOSYS
# define BIOLOGY
# define ANA_BIOLOGY
# define CARBON_ISOTOPE  /*Original CPP flag */
/*# define NUTRIENTS*/   /*Original CPP flag */
/*# define DENITRIFICATION*/
/*# define BIO_SEDIMENT*/
# define DIAGNOSTICS_BIO
#endif

