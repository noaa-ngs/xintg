#ifndef BILIN_H
#define BILIN_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: bilin.h 35300 2010-06-11 13:17:09Z Srinivas.Reddy $	20$Date: 2009/05/15 15:11:23 $ NGS"

// ----- standard library --------------------------------------------
#include <stdio.h>

// ----- classes, structs, types -------------------------------------
#include "grid_header.h"

// ----- functions ---------------------------------------------------
double bilin(double latdd, double londd,
             FILE* vec_ifp[50], GRID_HEADER vec_hdr[50], int kk);

#endif //~BILIN_H

