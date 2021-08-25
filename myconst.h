/* SccsID = "$Id: myconst.h 62154 2012-04-25 20:11:07Z Srinivas.Reddy $	%P%"  */
#ifndef __MYCONST__
#define __MYCONST__
#pragma ident "$Id: myconst.h 62154 2012-04-25 20:11:07Z Srinivas.Reddy $      20$Date: 2010/06/03 15:23:32 $ NGS"

#include <stdio.h>
#include <math.h>


#define MAX_MOSAIC                1100
#define MAX_MOSAIC_LIST           1100
#define MAX_MOSAIC_COUNTY        13000 // Max stations per County
#define MAX_MOSAIC_PROJ           1100 // Max stations per Project
#define MAX_MOSAIC_QUAD           1100 // Max stations per QUAD
#define MAX_MOSAIC_PROJ_INIT      2000 // Max stations per Project Initially
#define MAX_UNADJ_PROJ            5000
#define MAX_RADIUS_ITER             20
#define L_LENGTH                   255
#define MAX_CHAR                    80
#define MAX_ARC_CHAR              1024
#define MAX_CMD                    255
#define MAX_BLOCK                  300
#define MAX_ANT_BLOCK             1500
#define MAX_PH_PAT_BLOCK           400
#define MAX_EDIT_BLOCK              15
#define MAX_NET_STATUS              20
#define MAX_HISTORY                100
#define MAX_DESC                  1000
#define MAX_DIR                   1000
#define MAX_DIST                   100
#define MAX_LEV                    100
#define MAX_VERT                   100
#define MAX_LIST                  4500 // Big for 17137
#define MAX_ASTRO                  300
#define MAX_GRP                   10000 // was 20000 for GPS1751/B
#define MAX_VECTORS               18000
#define MAX_CORR                  2000
#define MAX_COV                    600
#define MAX_CORR_LOAD            80000
#define MAX_COV_LOAD             12000
#define MAX_DEFL                    20
#define MAX_ELLIP                   20
#define MAX_HGZ                    100
#define MAX_BOXSCORE               100
#define MAX_STATES                 300
#define MAX_REGIONS               5000 // Must match MAX_STATES_REGIONS for pos_to_state.c
#define MAX_CENSUS                2000
#define MAX_ERRORS                 100
#define MAX_UIDS                  1200
#define MAX_MARK                 20000
#define MAX_ABBR                   600
#define MAX_ALIAS                   10
#define MAX_PARENT                  50
#define MAX_DOPPLER                 50
#define MAX_DISTRICT               200
#define MAX_BORDER                1000
#define MAX_QUAKE                   50
#define MAX_AIRPORT                100 // 100 dates per AID
#define MAX_VERTCON_GRIDS            3
#define MAX_HPGN_GRIDS              55
#define MXAREA                      55
#define MAX_CAR97_GRIDS              2
#define MAX_MEX97_GRIDS              3
#define MAX_EGM96_GRIDS              1
#define MAX_DCAR97_GRIDS             6
#define MAX_GEOID99_GRIDS           14
#define MAX_GEOID03_GRIDS           14
#define MAX_GEOID06_GRIDS            4
#define MAX_GEOID09_GRIDS           16
#define MAX_DEFLEC_GRIDS            14
#define MAX_NADCON_GRIDS            11
#define MAX_USERS                  200
#define MAX_COUNTIES_REGIONS      1000
#define MAX_STATES_REGIONS        5000
#define MAX_REGION_PER_ID           30
#define MAX_COUNTIES_PER_STATE    1000
#define MAX_VM                      10
#define IMAGES_PER_UID              10
#define PI                         (4.0*atan(1.0))
#define TO_RADIANS(x)             ((x)*PI/180.0)
#define FR_RADIANS(x)             ((x)*180.0/PI)
#define TO_DEC_DEG(deg,min,sec)   ((double)deg + (double)min/60.0 + sec/3600.0)
#define MILES_TO_METERS(x)        ((x)*1609.3472)
#define METERS_TO_MILES(x)        ((x)/1609.3472)
#define MM_TO_METERS(x)           ((x)/1000.0)
#define CM_TO_METERS(x)           ((x)/100.0)
#define sFT_TO_METERS(x)          ((x)/39.37*12.00)
#define INCH_TO_METERS(x)          ((x)/39.37/12.00*12.00)
#define E2                        .0066943800229034156
#define AX                         6378137.0
#define NAD83_AX                   6378137.0
#define NAD27_AX                   6378206.4;
#define NAD83_RF                   298.257222101;
#define NAD27_RF                   294.9786982138;
#define QUAL_TEXT_LEN              200

typedef struct {
   char adj_id[14];
   char adj_date[9];
   char adj_agency[7];
   char adj_type[2];
   char load_date[9];
 } ADJUSTMENTS;

#endif
