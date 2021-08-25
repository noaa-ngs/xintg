#ifndef INTERG_IDW_H
#define INTERG_IDW_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: interg_idw.h 35302 2010-06-11 13:17:17Z Srinivas.Reddy $	20$Date: 2009/05/26 18:01:37 $ NGS"

// ----- standard library --------------------------------------------
#include <stdio.h>

// ----- classes, structs, types -------------------------------------
#include "grid_header.h"

// ----- functions ---------------------------------------------------
double interg_idw(double xlat, double xlon,
              GRID_HEADER  vec_hdr[50], FILE* vec_ifp[50], int kk);

#endif //~INTERG_IDW_H

