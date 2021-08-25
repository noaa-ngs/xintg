// %P%
// ----- constants ---------------------------------------------------
static const char SCCSID[]="$Id: ff2.c 43803 2010-07-12 16:43:46Z bruce.tran $	20$Date: 2009/05/27 15:01:21 $ NGS";
static const int  DEBUG = 0;           // diagnostics print if != 0

// ----- standard library --------------------------------------------
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>   // fabs()

// ----- classes, structs, types -------------------------------------
// ----- functions ---------------------------------------------------
#include "ff2.h"
#include "trim_c.h"


int ff2(char* card, double* lat_dd, double* lon_dd, char* text) {
/*******************************************************************************
* Reads a "Free Format, Type 2" record and 
* and parses latitude, longitude, descriptive text.
*   in - card : "Free Format, Type 2" record: 1-32 = lat, lon; 33-72 = text
*   in - lat_dd : latitude,         empty
*   in - lon_dd : longitude,        empty
*   in - text   : descriptive test, empty
*   out- lat_dd : latitude,         filled, decimal degrees
*   out- lon_dd : longitude,        filled, decimal degrees
*   out- text   : descriptive test  filled, char vector  (may be empty)
*   ret- zero for success
*******************************************************************************/
    char   SEPCHARS[] = " ";    // card is space delimited
    char   card_nums[42];
    char   card_text[42];
    char   numbers[7][30];
    char*  token;
    double latd=0.0, latm=0.0, lats=0.0;
    double lond=0.0, lonm=0.0, lons=0.0;
    double sign=1.0;
    int    ii = 0;

    // Initialize local variables
    for (ii = 0; ii < 7; ++ii) {
        strncpy(numbers[ii], "\0", 30);
    }

    strncpy(card_nums, "\0", 42);
    strncpy(text,      "\0", 42);

    strncpy(card_nums, &card[0], 32);
    strcpy( text,      &card[32]);

    trim_c(card_nums, 'b');
    // trim_c(text,      'b');


    // -----------------------------------------------------
    // Tokenize char* card_nums into char* numbers[]
    // -----------------------------------------------------
    ii = 0;
    if (strlen(card_nums) > 0) {
        token = strtok(card_nums, SEPCHARS);
        if (DEBUG > 0) { printf("token(1) = %s\n", token); }
        while (token != NULL) {
            strcpy(numbers[ii], token);
            if (DEBUG > 0) { printf("numbers[%d] = %s\n", ii, numbers[ii]); }
            ++ii;
            token = strtok(NULL, SEPCHARS);
        }
    }

    // Test for negative (South) latitude
    if (numbers[0][0] == '-') {
        sign = -1.0;
        numbers[0][0] = ' ';
    }

    switch (ii) {
    case 2:
        latd = atof( numbers[0] );
        lond = atof( numbers[1] );

        if (latd < 0) {           // Test for negative (South) latitude
            sign  = -1.0;
            latd *= -1.0;
        }
        *lat_dd = latd * sign;
        *lon_dd = lond;
        break;

    case 4:
        latd = atof( numbers[0] );
        latm = atof( numbers[1] );
        lond = atof( numbers[2] );
        lonm = atof( numbers[3] );

        if (latd < 0) {           // Test for negative (South) latitude
            sign  = -1.0;
            latd *= -1.0;
        }
        *lat_dd = (latm/60.0 + latd) * sign;
        *lon_dd = (lonm/60.0 + lond);
        break;

    case 6:
        latd = atof( numbers[0] );
        latm = atof( numbers[1] );
        lats = atof( numbers[2] );
        lond = atof( numbers[3] );
        lonm = atof( numbers[4] );
        lons = atof( numbers[5] );

        if (latd < 0) {           // Test for negative (South) latitude
            sign  = -1.0;
            latd *= -1.0;
        }
        *lat_dd = ((lats/60.0 + latm)/60.0 + latd) * sign;
        *lon_dd = ((lons/60.0 + lonm)/60.0 + lond);
        break;

    default:
        *lat_dd = -999.0;
        *lon_dd = -999.0;
        break;

    }//~switch

    return(0);

}//~ff2

