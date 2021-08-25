// %P%
// ----- constants ---------------------------------------------------
static const char SCCSID[]="$Id: intro.c 43798 2010-07-12 16:43:28Z bruce.tran $	20$Date: 2009/05/15 15:18:41 $ NGS";

// ----- standard library --------------------------------------------
#include <stdio.h>

// ----- classes, structs, types -------------------------------------
// ----- functions ---------------------------------------------------
#include "intro.h"


void intro(const char* PGMVER, const char* PGMDAT) {
/*******************************************************************************
* Print out introductory screen and disclaimer
* 
*******************************************************************************/

    printf("\n\
                             Welcome to the \n\
                       National Geodetic Survey's \n\
                               INTG PROGRAM \n\
                      (INTerpolate from Geoid files). \n\
            \n\
                                Version: %s \n\
                       Date(ccyy/mm/dd): %s \n\
                         (Hit RETURN to continue) \n", PGMVER, PGMDAT);

    getchar();

    printf("\n\
 ----------------------------------------------------------------------\n\
\n\
                             DISCLAIMER\n\
\n\
 This program and supporting information is furnished by the government of\n\
 the United States of America, and is accepted/used by the recipient with\n\
 the understanding that the U. S. government makes no warranties, express or\n\
 implied, concerning the accuracy, completeness, reliability, or suitability\n\
 of this program, of its constituent parts, or of any supporting data.\n\
\n\
 The government of the United States of America shall be under no liability\n\
 whatsoever resulting from any use of this program.  This program should\n\
 not be relied upon as the sole basis for solving a problem whose incorrect\n\
 solution could result in injury to person or property.\n\
\n\
 This program is the property of the government of the United States of\n\
 America. Therefore, the recipient further agrees not to assert proprietary\n\
 rights therein and not to represent this program to anyone as being other\n\
 than a government program.\n\
\n\
               (Hit RETURN to continue)\n");

    getchar();

    printf("\n\
 ----------------------------------------------------------------------\n\
\n\
 The flag for bad data is '-999' which will mostly appear if the input \n\
 latitude/longitude of a point falls outside of the borders of available files.\n\
 Note the longitude can be entered in either positive EAST or\n\
 as positive WEST.\n\
\n\
 Allowable formats for latitude and longitude, both keyboard and by file are:\n\
  1) Integer degrees                  2) Decimal degrees\n\
  3) Int Deg & Integer Minutes        4) Int Deg & Decimal Minutes\n\
  5) Int Deg, Int Min, & Integer Sec  6) Int Deg, Int Min & Decimal Seconds\n\
 Space delimited values only.  No commas.\n\
 Negative values are allowed, but only with the negative  sign preceding\n\
 and touching degrees (i.e. -100 42 37.3)\n\
\n\
              (Hit RETURN to continue)\n");

    getchar();

}//~intro

