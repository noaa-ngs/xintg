// %P%
// ----- constants ---------------------------------------------------
static const char SCCSID[]="$Id: run_bbk.c 64749 2012-08-27 15:55:58Z Srinivas.Reddy $	20$Date: 2010/06/21 09:13:44 $ NGS";
static const int  DEBUG = 0;           // diagnostics print if != 0

// ----- standard library --------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ----- classes, structs, types -------------------------------------
#include "grid_header.h"

// ----- functions ---------------------------------------------------
#include "run_bbk.h"
#include "getheaders.h"
#include "bb80ll.h"
#include "interg.h"
#include "which1.h"


void run_bbk(FILE* ifp, FILE* ofp, FILE** vec_ifp, GRID_HEADER vec_hdr[50], 
             char vec_fnames[50][256], int nfiles, int imodel) {
/*******************************************************************************
*   in - ifp    : input bluebook file
*   in - ofp    : output bluebook file, optional
*   in - vecIos : vector of file ptrs, pointing to grid data files
*   in - imodel : user specified geoid model
*  INPUT FILE, BLUE BOOK FORMAT
* Spin through the blue book file, stopping at each * *80* record, 
* and while there, look forward to the next records to see 
* if there is an *86* * record.  
* If so,  overwrite the geoid field in it with an interpolated value.
* If not, create an *86* record.
*
* Read one line, and if it is NOT an *84* record, 
* replicate it immediately into the output file (if we're using one).
* If the record is not an *80* record, 
* come back and get a new record.
*******************************************************************************/
    const char b40[]="                                        ";
    int jjj   = 0;
    int kount = 0;
    int keep  = 0;
    double xn =   -90.0;
    double xs =    90.0;
    double xw =   360.0;
    double xe =     0.0;
    double min=  1000.0;
    double max= -1000.0;
    double ave= 0.0;
    double std= 0.0;
    double rms= 0.0;
    double minlat=   90.0;
    double maxlat=  -90.0;
    double minlon=  360.0;
    double maxlon=    0.0;
    double fact;
    char   card[90];
    char   card2[90];
    char   card86[90];
    double xlat;           // lat of point in question - decdeg
    double xlon;           // lon of point in question - decdeg
    double val;            // interpolated value
    int    ight;           // geoid height * 1000. for output
    char   ight_c[12];     // geoid height * 1000. for output, string
    // ***** int    nfiles = 0;
    int    kk     = 0;

    // getheaders( vec_ifp, vec_hdr );

    strncpy(ight_c, "\0", 12);
    strncpy(card,   "\0", 90);
    strncpy(card2,  "\0", 90);
    strncpy(card86, "\0", 90);

    // Spin through the blue book file, stopping at each *80* record,
    // and while there, look forward to the next records 
    // to see if there is an *86* record. 
    // If so, overwrite the geoid field in it with an interpolated value.
    // If not, create an *86* record.

    while( fgets(card, 90, ifp) ) {

        // Read one line, and if it is NOT an *84* record,
        // replicate it immediately into the output file (if we're using one).
        // If the record is not an *80* record, come back and get a new record.

        if( (strncmp(&card[7], "83", 2) == 0)
        ||  (strncmp(&card[7], "84", 2) == 0)) {
            continue;
        } else 
        if (ofp) {
            fprintf(ofp, "%s", card);
        }

        if (strncmp(&card[7], "80", 2) != 0) {
            continue;
        } else {
            // arriving here, we have an *80* record as 'card'
            bb80ll(card, &xlat, &xlon);

            if (DEBUG != 0) {
                printf("card = \n%s\n", card);
                printf("xlat = %lf  xlon = %lf\n", xlat, xlon);
                printf("card[55] = %c \n", card[55]);
                printf("card[68] = %c \n", card[68]);
            }

            // If the xlat/xlon values are bad, don't do an interpolation.
            // Just move on to the next record...any existing *86*
            // record associated with this erroneous *80* record will
            // remain unmodified.
            if (xlat == -999.0 || xlon == -999.0)  continue;

            if (card[55] == 'S' || card[55] == 's') xlat *= -1.0;

            if (card[68] == 'W' || card[68] == 'w') xlon = 360.0 - xlon;

            if (DEBUG != 0) {
                printf("xlat = %lf  xlon = %lf\n", xlat, xlon);
            }

            // Arriving here, the *80* record has a good lat/lon value
            ++kount;

            // The associated *86* record, if it exists, must
            // come AFTER the *80* record we're currently looking at, 
            // and the ONLY thing allowed between the 80 and 86 record 
            // are 81, 82 and 85 records (83 and 84 records will be deleted).

            jjj = 0;
            while (jjj == 0) {
                strncpy(card2, "\0", 90);
                fgets(card2, 90, ifp);

                if (DEBUG != 0) {
                    printf("card2 = \n%s\n", card2);
                }

                // POSSIBILITY 1: We find an *83* or *84* record
                // ACTION: Do not allow it to go into the output file, and
                //         go back for another "card2" value

                if( (strncmp(&card2[7], "83", 2) == 0) 
                ||  (strncmp(&card2[7], "84", 2) == 0) ) {
                    continue;
                }

                // POSSIBILITY 2: We find an *81*, *82* or *85* record
                // ACTION: Replicate this record into the output file
                //         and go back for another "card2" value

                if( (strncmp(&card2[7], "81", 2) == 0) 
                ||  (strncmp(&card2[7], "82", 2) == 0) 
                ||  (strncmp(&card2[7], "85", 2) == 0) ) {
                    fprintf(ofp, "%s", card2);
                    continue;
                }

                // POSSIBILITY 3: We find an *86* record
                // ACTION: a) Modify it if it's the right *86* record
                //         b) If it is the wrong *86* record, 
                //            delete it and make a NEW *86* record, and
                //            then go back and get another "card2" value

                if( (strncmp(&card2[7], "86", 2)) == 0)  {

                    if( (strncmp(&card[10], &card2[10], 4)) == 0)  {

                        if (DEBUG != 0) { printf("match 80 and 86"); }

                        // Arriving here, we've found the *86* record 
                        // associated with the *80* record

                        // Find which geoid file to use, based on lat/lon 
                        kk = which1(xlat, xlon, nfiles, kk, imodel, 
                               vec_fnames, vec_hdr, vec_ifp );
                        if (DEBUG != 0) {
                            printf("After which1 kk = %d\n", kk);
                        }

                        if (kk == -1) {
                            // If not in any grid area, set to -999
                            val = -999;
                        } else {
                            // Otherwise, do the interpolation

                            val = interg(xlat, xlon, vec_hdr, vec_ifp, kk);

                            keep +=   1;
                            ave  += val;
                            rms  += (val*val);
                            if (xlat > xn)  xn = xlat;
                            if (xlat < xs)  xs = xlat;
                            if (xlon > xe)  xe = xlon;
                            if (xlon < xw)  xw = xlon;
                            if (val < min)  {
                                min = val;
                                minlat = xlat;
                                minlon = xlon;
                            }
                            if (val > max) {
                                max = val;
                                maxlat = xlat;
                                maxlon = xlon;
                            }

                            ight = (int)floor((val * 1000.0) + 0.5);
                            sprintf(ight_c, "%7d", ight);

                            strncpy( &card2[35], ight_c, 7);

                            // Modify the  *86* record and write it out
                            // Find the right geoid code
                            if (imodel == 1)  card2[42] = 'T';   // GEOID99
                            if (imodel == 2)  card2[42] = 'U';   // GSSS99
                            if (imodel == 3)  card2[42] = 'W';   // GEOID03
                            if (imodel == 4)  card2[42] = 'X';   // USGG2003
                            if (imodel == 5)  card2[42] = 'Y';   // GEOID06
           // Not Used    if (imodel == *)  card2[42] = 'Z';   // USGG2006
                            if (imodel == 6)  card2[42] = '1';   // USGG2009
                            if (imodel == 7)  card2[42] = '2';   // GEOID09
                            if (imodel == 11) card2[42] = '4';   // USGG2012
                            if (imodel == 12) card2[42] = '5';   // GEOID12A

                            if (ofp) {
                                fprintf(ofp, "%s", card2);
                            }

                            // Get the next good *80* record in "card"
                            //    goto 137
                            ++jjj;
                            strncpy(card,  "\0", 90);
                            strncpy(card2, "\0", 90);
                            continue;
                        }

                    }//~if( strncmp(card2[11], card[11], 4) == 0 )  

                    // If this is NOT the associated *86* record, delete
                    // it and make a new *86* record
                    else {
                        if (DEBUG != 0) { printf("NOT match 80 and 86"); }

                        // Find geoid file to use, based on input lat/lon 
                        kk = which1(xlat, xlon, nfiles, kk, imodel, 
                               vec_fnames, vec_hdr, vec_ifp );

                        if (kk == -1) {
                            // If not in any grid area, set = -999
                            val = -999;
                        } else {
                            // Otherwise, do the interpolation

                            val = interg(xlat, xlon, vec_hdr, vec_ifp, kk);

                            keep +=   1;
                            ave  += val;
                            rms  += (val*val);
                            if (xlat > xn)  xn = xlat;
                            if (xlat < xs)  xs = xlat;
                            if (xlon > xe)  xe = xlon;
                            if (xlon < xw)  xw = xlon;
                            if (val < min)  {
                                min = val;
                                minlat = xlat;
                                minlon = xlon;
                            }
                            if (val > max) {
                                max = val;
                                maxlat = xlat;
                                maxlon = xlon;
                            }

                            // Make a NEW *86* record and write it out
                            strncpy(card86, "\0",     90);
                            strncpy(card86, card,      6);
                            strncat(card86, "*86*",    4);
                            strncat(card86, &card[10], 4);
                            strcat( card86, "                   ");
                            strcat( card86, "                   ");

                            ight = (int)floor((val * 1000.0) + 0.5);
                            sprintf(ight_c, "%7d", ight);

                            strncpy( &card86[35], ight_c, 7);

                            // Find the right geoid code
                            if (imodel == 1)  card86[42] = 'T';  // GEOID99
                            if (imodel == 2)  card86[42] = 'U';  // GSSS99
                            if (imodel == 3)  card86[42] = 'W';  // GEOID03
                            if (imodel == 4)  card86[42] = 'X';  // USGG2003
                            if (imodel == 5)  card86[42] = 'Y';  // GEOID06
           // Not Used    if (imodel == *)  card2[42] = 'Z';   // USGG2006
                            if (imodel == 6)  card86[42] = '1';  // USGG2009
                            if (imodel == 7)  card86[42] = '2';  // GEOID09
                            if (imodel == 11) card86[42] = '4';   // USGG2012
                            if (imodel == 12) card86[42] = '5';   // GEOID12A

                            if (ofp) {
                                fprintf(ofp, "%s\n", card86);
                            }

                            // Go get the next good *80* record in "card"
                            //     goto 137
                            ++jjj;
                            continue;

                        }//~if (kk == -1) 

                    }//~if (card2.compare(11,4, card,11,4) == 0)

                strncpy(card2, "\0", 90);
                }//~if (card2.compare(7, 2, "86") == 0)


                // POSSIBILITY 4: We find something besides an 81,82,
                //                83, 84, 85 or 86 record
                // ACTION: Create a new *86* based on the previous *80* record,
                //         and then transfer the new 'card2' into 'card' and
                //         go back just after the 'read' statement for 'card'

                if( (strncmp(&card2[7], "81", 2) != 0)
                &&  (strncmp(&card2[7], "82", 2) != 0)
                &&  (strncmp(&card2[7], "83", 2) != 0)
                &&  (strncmp(&card2[7], "84", 2) != 0)
                &&  (strncmp(&card2[7], "85", 2) != 0)
                &&  (strncmp(&card2[7], "86", 2) != 0) ) { 

                    // Find geoid file to use, based on input lat/lon 
                    kk = which1(xlat, xlon, nfiles, kk, imodel, vec_fnames, vec_hdr,vec_ifp);

                    if (kk == -1) {
                        // If not in any of our grid areas, set to -999
                        val = -999;
                    } else {
                        // Otherwise, do the interpolation

                        val = interg(xlat, xlon, vec_hdr, vec_ifp, kk);

                        keep += 1;
                        ave  += val;
                        rms  += val*val;
                        if (xlat > xn)  xn = xlat;
                        if (xlat < xs)  xs = xlat;
                        if (xlon > xe)  xe = xlon;
                        if (xlon < xw)  xw = xlon;
                        if (val < min) {
                            min = val;
                            minlat = xlat;
                            minlon = xlon;
                        }
                        if (val > max) {
                            max = val;
                            maxlat = xlat;
                            maxlon = xlon;
                        }
                    }

                    // Make a NEW *86* record
                    strncpy(card86, "\0",     90);
                    strncpy(card86, card,      6);
                    strncat(card86, "*86*",    4);
                    strncat(card86, &card[10], 4);
                    strcat( card86, b40);

                    // Find the right geoid code
                    if (imodel == 1)  card86[42] = 'T';  // GEOID99
                    if (imodel == 2)  card86[42] = 'U';  // GSSS99
                    if (imodel == 3)  card86[42] = 'W';  // GEOID03
                    if (imodel == 4)  card86[42] = 'X';  // USGG2003
                    if (imodel == 5)  card86[42] = 'Y';  // GEOID06
      // Not Used  if (imodel == *)  card2[42] = 'Z';   // USGG2006
                    if (imodel == 6)  card2[42] = '1';   // USGG2009
                    if (imodel == 7)  card2[42] = '2';   // GEOID09

                    ight = (int)floor((val * 1000.0) + 0.5);
                    sprintf(ight_c, "%7d", ight);

                    strncpy( &card86[35], ight_c, 7);

                    if (ofp) {
                       fprintf(ofp, "%s\n", card86);
                    }

                    // Now put the latest *80* record, which was in 'card2'
                    // into 'card' and go back to the top to search

                    strcpy(card, card2);

                }//~if ((card2.substr(7,2) != "81")


                // POSSIBILITY 5: We find the EOF
                // ACTION: Create a new *86* based on the previous *80* record
                //         and then go to the 'end report' phase

                // Find geoid file to use, based on input lat/lon 
                kk = which1(xlat, xlon, nfiles, kk, imodel, vec_fnames, vec_hdr,vec_ifp);

                // If the point isn't in any of our grid areas, set to -999
                if (kk == -1) {
                    val = -999;
                } else {
                    // Otherwise, do the interpolation

                    val = interg(xlat, xlon, vec_hdr, vec_ifp, kk);

                    keep +=   1;
                    ave  += val;
                    rms  += (val*val);
                    if (xlat > xn)  xn = xlat;
                    if (xlat < xs)  xs = xlat;
                    if (xlon > xe)  xe = xlon;
                    if (xlon < xw)  xw = xlon;
                    if (val < min) {
                        min = val;
                        minlat = xlat;
                        minlon = xlon;
                    }
                    if (val > max) {
                        max = val;
                        maxlat = xlat;
                        maxlon = xlon;
                    }
                }

                // Make a NEW *86* record
                strncpy(card86, "\0",     90);
                strncpy(card86, card,      6);
                strncat(card86, "*86*",    4);
                strncat(card86, &card[10], 4);
                strcat( card86, b40);

                // Find the right geoid code
                if (imodel == 1)  card2[42] = 'T';     // GEOID99
                if (imodel == 2)  card2[42] = 'U';     // GSSS99
                if (imodel == 3)  card2[42] = 'W';     // GEOID03
                if (imodel == 4)  card2[42] = 'X';     // USGG2003
                if (imodel == 5)  card2[42] = 'Y';     // GEOID06
  // Not Used  if (imodel == *)  card2[42] = 'Z';   // USGG2006
                if (imodel == 6)  card2[42] = '1';   // USGG2009
                if (imodel == 7)  card2[42] = '2';   // GEOID09

                ight = (int)floor((val * 1000.0) + 0.5);
                sprintf(ight_c, "%7d", ight);

                strncpy( &card86[35], ight_c, 7);

                if (ofp) {
                   //if (strncmp(ight_c, "-999000",7)!=0) {
                      fprintf(ofp, "%s\n", card86);
                   //}
                }

                strncpy(card,  "\0", 90);
                strncpy(card2, "\0", 90);
                jjj=1;
            }//~while (jjj == 0)

        }//~if ( card.compare(7, 2, "80") != 0)

    }//~while( fgets(card, 80, ifp) ) {

    // When finished, give a little report to the screen and end program.
    ave /= keep;
    rms = sqrt(rms/keep);
    if (keep > 1) {
        fact = ((double)(keep))/((double)(keep-1));
        std = sqrt( fact * ((rms*rms) - (ave*ave)) );
    } else {
        std = 0;
    }

    printf(\
" ----------------------------------------------------------------------\n\
 FINAL REPORT: \n\
  Number of points input: %d \n\
  Number of good points : %d \n\
  Northernmost Latitude : %10.6lf \n\
  Southernmost Latitude : %10.6lf \n\
  Westernmost Longitude : %10.6lf \n\
  Easternmost Longitude : %10.6lf \n", kount, keep, xn, xs, xw, xe);

    printf(\
"  Minimum Geoid Height  : %8.3lf \n\
    Lat/Lon of Minimum  : %10.6lf  %10.6lf \n\
  Maximum Geoid Height  : %8.3lf \n\
    Lat/Lon of Maximum  : %10.6lf  %10.6lf \n\
  Average Geoid Height  : %8.3lf \n\
  Standard Deviation    : %8.3lf \n\
  Root Mean Square      : %8.3lf \n", 
    min, minlat, minlon, max, maxlat, maxlon, ave, std, rms);

    return;

}//~run_bbk

