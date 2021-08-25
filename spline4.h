#ifndef SPLINE4_H
#define SPLINE4_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: spline4.h 35297 2010-06-11 13:16:57Z Srinivas.Reddy $	20$Date: 2009/05/15 14:00:54 $ NGS"

// ----- standard library --------------------------------------------
#include <stdio.h>

// ----- classes, structs, types -------------------------------------
#include "grid_header.h"

// ----- functions ---------------------------------------------------
double spline4(double latdd,      double londd,
               FILE* vec_ifp[50], GRID_HEADER vec_hdr[50], int kk);

#endif //~SPLINE4_H

