// %P%
// ----- constants ---------------------------------------------------
static const char SCCSID[]="$Id: bb80ll.c 43793 2010-07-12 16:43:10Z bruce.tran $	20$Date: 2009/05/15 15:10:14 $ NGS";

// ----- standard library --------------------------------------------
#include <string.h>

// ----- classes, structs, types -------------------------------------
#include "dms_dd.h"


int  bb80ll( char* card, double* lat, double* lon ) {
/*******************************************************************************
* Parse *80 record, returns lat and lon by reference
*   in - card : rec_in
*   out- lat  : North latitude,   decimal degrees
*   out- lon  : West  longitude,  decimal degrees
* 
*******************************************************************************/
    char latc[30];
    char lonc[30];

    strncpy( latc, "\0", 30);
    strncpy( lonc, "\0", 30);

    strncpy( latc, &card[55], 1);           // H            at col 56
    strncat( latc, &card[44], 6);           // DDMMSS       at col 45
    strcat(  latc, ".");
    strncat( latc, &card[50], 5);           // sssss        at col 51

    strncpy( lonc, &card[68], 1);           // H            at col 69
    strncat( lonc, &card[56], 7);           // DDDMMSS      at col 57
    strcat(  lonc, ".");
    strncat( lonc, &card[63], 5);           // sssss        at col 64

    *lat = dms_dd( latc );
    *lon = dms_dd( lonc );

    return( 0 );

}//~bb80ll

