#ifndef IDW_INTERP_H
#define IDW_INTERP_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: idw_interp.h 35297 2010-06-11 13:16:57Z Srinivas.Reddy $	20$Date: 2009/05/15 14:00:54 $ NGS"

// ----- standard library --------------------------------------------
#include <stdio.h>

// ----- classes, structs, types -------------------------------------
#include "idw_point.h"

// ----- functions ---------------------------------------------------
double idw_interp(double latdd, double londd, IDW_POINT *points, int num);

#endif //~IDW_INTERP_H

