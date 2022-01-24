// %P%
// ----- constants ---------------------------------------------------
static const char SCCSID[]="$Id: c2v.c 35388 2010-06-11 14:02:16Z Srinivas.Reddy $	20$Date: 2009/05/15 13:57:26 $ NGS";
static const int  DEBUG = 0;           // diagnostics print if != 0

// ----- standard library --------------------------------------------
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>    // fabs

// ----- classes, structures, types ----------------------------------
// ----- functions ---------------------------------------------------
#include "c2v.h"
//#include "trim_c.h"


double c2v(char* card, int ilatlon) {
/*******************************************************************************
* Function to convert a 40 character char* into a value.
*
* Reads a string, tokenizes, forms and returns a value of type double (decdeg)
* The string must consist of 1,2 or 3 numbers, meaning
* Degrees, Degrees&Minutes or Deg, Min & Sec
*
* The allowable characters in "card" are spaces, numbers 1 through 9,
* decimal points, and one leading negative sign.
*
* Returned values of "value" will be forced into:
*   ilatlon = 1 =>     -90 < value < +90
*   ilatlon = 2 =>       0 < value < +360
*
* This has nothing to do with Positive Longitude conventions.
* It is purely a character-to-value conversion.
* Longitude conventions are handled outside this function.
*******************************************************************************/
    char   SEPCHARS[] = " ";    // card is space delimited
    char   numbers[4][30];
    char*  token;
    double value=0.0;   // return value
    double deg=0.0;
    double min=0.0;
    double sec=0.0;
    double sign = 1.0;
    int    ii = 0;

    if (DEBUG > 0) { printf("card = %s\n", card); }

    for (ii = 0; ii < 4; ++ii) {
        strncpy(numbers[ii], "\0", 30);
    }

    trim_c(card, 'B');

    // Tokenize char* card into char* vector numbers
    ii = 0;
    if (strlen(card) > 0) {
        token = strtok(card, SEPCHARS);
        if (DEBUG > 0) { printf("token(1)  = %s\n", token); }
        while (token != NULL) {
            strcpy(numbers[ii], token);
            if (DEBUG > 0) { printf("numbers[%d] = %s\n", ii, numbers[ii]); }
            ++ii;
            token = strtok(NULL, SEPCHARS);
        }
    } else {
        printf("Error: Empty string into function c2v\n");
    }

    switch (ii) {
    case 1:
        value = atof( numbers[0] );
        break;

    case 2:
        deg = atof( numbers[0] );
        min = atof( numbers[1] );

        if (deg < 0) {
            sign = -1.0;
            deg = fabs(deg);
        }

        value  = min/60.0 + deg;
        value *= sign;
        break;

    case 3:
        deg = atof( numbers[0] );
        min = atof( numbers[1] );
        sec = atof( numbers[2] );

        if (deg < 0) {
            sign = -1.0;
            deg = fabs(deg);
        }

        value  = (sec/60.0 + min)/60.0 + deg;
        value *= sign;
        break;

    default:
        value = -999.0;

    }//~switch

    // Normalize the lat or long value
    if (fabs(value) < 998.0) {
        if (ilatlon == 1) {   // is lat
            if (value < -90.0) value = -999.0;
            if (value >  90.0) value = -999.0;
        } else
        if (ilatlon == 2) {   // is lon
            while (value <   0.0)  value += 360.0;
            while (value > 360.0)  value -= 360.0;
        } else {
            printf("Error: in function c2v, ilatlon out of range = %d\n",
                ilatlon);
        }
    }

    return( value );

}//~c2v

