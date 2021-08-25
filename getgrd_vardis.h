#ifndef GETGRD_VARDIS_H
#define GETGRD_VARDIS_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: getgrd_vardis.h 35309 2010-06-11 13:17:45Z Srinivas.Reddy $	20$Date: 2010/03/24 15:22:51 $ NGS"

// ----- standard library --------------------------------------------
#include <stdio.h>

// ----- classes, structs, types -------------------------------------
// ----- functions ---------------------------------------------------
void getgrd_vardis(int  imodel, char* dirnam, int is_subr, int* nfiles,
                 int*  nvff,                 int*  ndff,  
                 char  var_fnames[50][256], 
                 FILE* var_ifp[50] );

#endif //~GETGRD_VARDIS_H

