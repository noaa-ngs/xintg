// %P%
// ----- constants ---------------------------------------------------
static const char SCCSID[]="$Id: ff2out.c 67515 2012-12-20 20:57:13Z stephen.gilbert $	20$Date: 2010/03/01 13:24:03 $ NGS";

// ----- standard library --------------------------------------------
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// ----- classes, structs, types -------------------------------------
#include "dataset1.h"

// ----- functions ---------------------------------------------------
#include "ff2out.h"
#include "dd_dms.h"


int ff2out(FILE* ofp, DATASET1 vec_data, double geoidHt, int imodel,
           double stddev, double distance) {
/*******************************************************************************
* Writes a "Free Format, (For Geoid) Type 2" record to an output file
*   in - ofp      : pointer to output file
*   in - vec_data : vector containing input data
*   in - geoidHt  : interpolated geoid height
*   out- nothing
*   ret- zero always
* was: string ff2out( iinput,outfil,lout,card, poseast,xlat,xlon,val ) {
*******************************************************************************/
    char  lat_dms[20];
    char  lon_dms[20];

    char  latdeg_c[12];
    char  latmin_c[12];
    char  latsec_c[12];
    char  londeg_c[12];
    char  lonmin_c[12];
    char  lonsec_c[12];
    char  cstddev[9];
    char  cdist[9];

    int   latdeg, latmin;
    int   londeg, lonmin;
    double latsec, lonsec;

    // Initialize local variables
    strncpy(lat_dms,  "\0", 20);
    strncpy(lon_dms,  "\0", 20);
    strncpy(latdeg_c, "\0", 12);
    strncpy(latmin_c, "\0", 12);
    strncpy(latsec_c, "\0", 12);
    strncpy(londeg_c, "\0", 12);
    strncpy(lonmin_c, "\0", 12);
    strncpy(lonsec_c, "\0", 12);
    strncpy(cstddev,  "\0",  9);
    strncpy(cdist,    "\0",  9);

    if (vec_data.lat == -999  ||  vec_data.lon == -999) {
        latdeg = 99;
        latmin = 99;
        latsec = 99.99999;

        londeg = 99;
        lonmin = 99;
        lonsec = 99.9999;
    } else {
        dd_dms( vec_data.lat, lat_dms );

        if (vec_data.poseast == 0) {     // 0 := false => is West lon)
            dd_dms( (360.0 - vec_data.lon), lon_dms );
        } else {
            dd_dms( vec_data.lon, lon_dms );
        }

        strncpy( latdeg_c, &lat_dms[0], 3);
        strncpy( latmin_c, &lat_dms[3], 2);
        strncpy( latsec_c, &lat_dms[5], 8);

        strncpy( londeg_c, &lon_dms[0], 3);
        strncpy( lonmin_c, &lon_dms[3], 2);
        strncpy( lonsec_c, &lon_dms[5], 8);

        latdeg = atoi( latdeg_c );
        latmin = atoi( latmin_c );
        latsec = atof( latsec_c );

        londeg = atoi( londeg_c );
        lonmin = atoi( lonmin_c );
        lonsec = atof( lonsec_c );

        if ((latsec + 0.000001) >= 60.0) {
            latsec -= 60.0;
            ++latmin;
        }
        if (latmin >= 60) {
            latmin -= 60;
            ++latdeg;
        }

        if ((lonsec + 0.000001) >= 60.0) {
            lonsec -= 60.0;
            ++lonmin;
        }
        if (lonmin >= 60) {
            lonmin -= 60;
            ++londeg;
        }
    }

    if ( stddev == -999. ) {
       strcpy(cstddev," UNAVAIL");
    }
    else {
       sprintf(cstddev,"%8.4lf",stddev);
    }

    if ( distance == -999. ) {
       strcpy(cdist," UNAVAIL");
    }
    else {
       sprintf(cdist,"%8.3lf",distance);
    }

    // -----------------------------------------------
    // Print to an output file
    // -----------------------------------------------
    if ( imodel == 12 ) 
       fprintf(ofp, "%3d %2d %8.5lf %3d %2d %8.5lf %8.3lf %s %s %s\n",
               latdeg, latmin, latsec, londeg, lonmin, lonsec,
               geoidHt, cstddev, cdist, vec_data.text);
    else
       fprintf(ofp, "%3d %2d %8.5lf %3d %2d %8.5lf %8.3lf%s\n",
               latdeg, latmin, latsec, londeg, lonmin, lonsec,
               geoidHt, vec_data.text);

    return( 0 );

}//~ff2out

