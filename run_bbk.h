#ifndef RUN_BBK_H
#define RUN_BBK_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: run_bbk.h 61167 2012-03-06 21:09:48Z bill.waickman $	20$Date: 2009/05/15 15:19:29 $ NGS"

// ----- standard library --------------------------------------------
#include <stdio.h>

// ----- classes, structs, types -------------------------------------
#include "grid_header.h"

// ----- functions ---------------------------------------------------
void run_bbk(FILE* ifp, FILE* ofp, FILE* vec_ifp[50], GRID_HEADER vec_hdr[50], 
             char vec_fnames[50][256], int nfiles, int imodel);

#endif //~RUN_BBK_H

