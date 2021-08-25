#ifndef FF2OUT_H
#define FF2OUT_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: ff2out.h 67515 2012-12-20 20:57:13Z stephen.gilbert $	20$Date: 2009/07/06 17:55:27 $ NGS"

// ----- standard library --------------------------------------------
#include <stdio.h>

// ----- classes, structs, types -------------------------------------
#include "dataset1.h"

// ----- functions ---------------------------------------------------
int ff2out(FILE* ofp, DATASET1 vec_data, double geoidHt, int imodel,
           double stddev, double distance);

#endif //~FF2OUT_H

