// %P%
// ----- constants ---------------------------------------------------
static const char SCCSID[]="$Id: expform.c 43801 2010-07-12 16:43:38Z bruce.tran $	20$Date: 2009/05/22 13:12:39 $ NGS";

// ----- standard library --------------------------------------------
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// ----- classes, structs, types -------------------------------------
// ----- functions ---------------------------------------------------
#include "expform.h"


void expform( void ) {
/*******************************************************************************
* Prints an explanation of the various formats to the display
* 
*******************************************************************************/
    char  keyb[11];
    char  cformx[40];
    int   iformx =0;
    int   iii    =0;

    strncpy(keyb,   "\0", 11);
    strncpy(cformx, "\0", 40);


    while (iii == 0) {
    // -----------------------------------------------
    // Find out which format needs explaining
    // -----------------------------------------------
    printf("\
Which format would you like explained? \n\
  1 = Free Format Type 1               \n\
  2 = Free Format Type 2               \n\
  3 = NGS Blue Book Format             \n\
 99 = Return to the main program       \n\
  0 = End Program                      \n\
  -> ");

    fgets(cformx, 40, stdin);
    iformx = atoi( cformx );

    switch (iformx) {
    case 0: 
        exit(0);
    case 99: 
        ++iii;
        break;

    case 1: 
    // -------------------------------------------
    // Explain Free Format (For Geoid) Type 1
    // -------------------------------------------
        printf("\
----------------------------------------------------------------------\n\
 1 = Free Format (For Geoid) Type 1 : \n\
This is an ASCII format, where the first 40 characters  \n\
of each record may contain the station name or be blank.\n\
The rest of the record (characters 41-80) must contain  \n\
the latitude and longitude in one of three formats:     \n\
  (a) decimal/integer  degrees  (two numbers)           \n\
  (b) integer degrees, decimal/integer minutes (four numbers) \n\
  (c) integer degrees, integer minutes, decimal/integer seconds \
(six numbers) \n\
The latitude must be positive North. \n\
The longitude sign convention used is specified by the users \
during the running of this program. \n\n\
                         (Hit RETURN to continue)\n");

        getchar();

        printf("\
The following records are examples of acceptable input: \n\n\
         1         2         3         4         5         6         7         8\n\
12345678901234567890123456789012345678901234567890123456789012345678901234567890\n\
                                        <------------ Columns 41-80 ----------->\n\
BM1027                                  50 285 \n\
XX8063                                          45.3 282 \n\
                 AIRPORT SC1               47 18.55 252 14 \n\
                                         30 17 270 59 \n\
           CORNER MARKER                 48 14 17 239 52 18.753 \n\
 DATA PT #4077                             42 42 42.4242 255 55 55.55555\n\n\
                         (Hit RETURN to continue)\n");

        getchar();
        break;

    case 2: 
    // -------------------------------------------
    // Explain Free Format (For Geoid) Type 2
    // -------------------------------------------
        printf("\
----------------------------------------------------------------------\n\
 2 = Free Format (For Geoid) Type 2 : \n\
This is an ASCII format, where the last 40 characters \n\
of each record may contain the station name or be blank.\n\
The first 32 characters of the record (characters 1-32)\n\
must contain the latitude and longitude in one of three \
 formats:\n\
   (a) decimal/integer  degrees  (two numbers)\n\
   (b) integer degrees, decimal/integer minutes  (four numbers)\n\
   (c) integer degrees, integer minutes, decimal/integer \
 seconds  (six numbers)\n\
The latitude must be positive North.\n\
The longitude sign convention used is specified by the users \
during the running of this program.\n\n\
                         (Hit RETURN to continue)\n");

        getchar();

        printf("\
 The following records are examples of acceptable input: \n\n\
         1         2         3         4         5         6         7         8\n\
12345678901234567890123456789012345678901234567890123456789012345678901234567890\n\
<-------- Columns 1-32 --------> \n\
50 285                          BM1027 \n\
         45.3 282               XXX8063 \n\
  47 18.55 252 14                               AIRPORT SC1 \n\
30 17 270 59 \n\
48 14 17 239 52 18.753                    CORNER MARKER \n\
  42 42 42.4242 255 55 55.55555 DATA PT #4077 \n\n\
                         (Hit RETURN to continue)\n");

        getchar();
        break;

    case 3: 
    // -------------------------------------------
    // Explain NGS Blue Book Format
    // -------------------------------------------
        printf("\
----------------------------------------------------------------------\n\
 3 = NGS Blue Book Format : \n\n\
 NGS Horizontal (Vol 1) Blue Book format - \
*80* (Control Point) and *86* (Combo Height) Records.\n\
INTG uses information from the *80* and *86* records in Blue \
Book files. \n\
The Station Serial Number(SSN), \
and the Latitude and the Longitude\n\
are read from the *80* records.  \n\
The rest of the record is ignored. \n\
Only the *86* record is modified by INTG - \
The *80* records and all other records are \
passed through without change to the output file.\n\
Any *84* records are deleted.\n\
Information in cols 36-43 in the *86* record is created by INTG.\n\
\n\
For more information on this format, please refer to: \n\
   \'Input Formats and Specifications of the \
National Geodetic Survey Data Base\' \n\
   \'Volume 1. Horizontal Control Data\'. \n\
Published by the Federal Geodetic Control Committee \
in September 1994 and \n\
available from: the National Geodetic Survey, NOAA, \
 Silver Spring, MD 20910. \n\n\
                         (Hit RETURN to continue)\n");

        getchar();

        printf("\
 The following input example is *80* and *86* records \
from a Blue Book file: \n\n\
004560*80*0096KNOXVILLE CT HSE              411906578  N0930548534  W 277  MIA33\n\
004565*86*0096 any geoid values in here     will be overwritten if ssn matches  \n\
                         (Hit RETURN to continue)\n");

        getchar();
        break;

    default:
        printf("Error: Invalid response\n");

    }//~switch 
    }//~while(iii)

}//~expform

