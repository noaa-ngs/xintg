#ifndef _GETGRD_GEOID_H_
#define _GETGRD_GEOID_H_

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: getgrd_geoid.h 35307 2010-06-11 13:17:37Z Srinivas.Reddy $	20$Date: 2010/02/22 09:41:14 $ NGS"

// ----- standard library --------------------------------------------
#include <stdio.h>

// ----- classes, structs, types -------------------------------------
// ----- functions ---------------------------------------------------
void getgrd_geoid(int imodel, char* dirnam, int is_subr, int* nfiles, int* nff,
                  char vec_fnames[50][256], FILE* vec_ifp[50]);

#endif //~_GETGRD_GEOID_H_

