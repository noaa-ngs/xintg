#ifndef SPLINE6_H
#define SPLINE6_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: spline6.h 35298 2010-06-11 13:17:01Z Srinivas.Reddy $	20$Date: 2009/05/15 14:01:00 $ NGS"

// ----- standard library --------------------------------------------
#include <stdio.h>

// ----- classes, structs, types -------------------------------------
#include "grid_header.h"

// ----- functions ---------------------------------------------------
double spline6(double latdd,      double londd,
               FILE* vec_ifp[50], GRID_HEADER vec_hdr[50], int kk);

#endif //~SPLINE6_H

