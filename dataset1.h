#ifndef DATASET1_H
#define DATASET1_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: dataset1.h 37479 2010-06-16 18:36:45Z bruce.tran $	20$Date: 2010/03/01 11:05:11 $ NGS"

// ----- standard library --------------------------------------------
// ----- classes, structs, types -------------------------------------
typedef struct {
    double lat;     // decimal degrees
    double lon;     // decimal degrees, pos east
    int    poseast; // Cstd: 0 = false(is West lon); !0 = not false(is East lon)
    char text[50];  // ff1: 1st 40 char is field  ff2: last 40 char is field
} DATASET1;

#endif //~DATASET1_H

