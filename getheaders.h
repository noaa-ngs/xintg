#ifndef GETHEADERS_H
#define GETHEADERS_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: getheaders.h 35294 2010-06-11 13:16:45Z Srinivas.Reddy $	20$Date: 2009/05/15 14:00:14 $ NGS"

// ----- standard library --------------------------------------------
#include "stdio.h"

// ----- classes, structs, types -------------------------------------
#include "grid_header.h"            // struct

// ----- functions ---------------------------------------------------
int getheaders( FILE** vec_ifp, GRID_HEADER* vec_hdr, int nfiles );

#endif //~GETHEADERS_H

