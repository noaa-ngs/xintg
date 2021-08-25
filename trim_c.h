#ifndef TRIM_C_H
#define TRIM_C_H

#pragma ident "$Id: trim_c.h 63417 2012-06-25 20:07:13Z Srinivas.Reddy $	20$Date: 2010/06/03 15:23:52 $ NGS"
#include <string.h>

//-----------------------------------------------------------------------------
// Removes leading and/or trailing white space chars from the input string.
// The second argument indicates whether to trim
// blanks from both ends (the default), the left side, or the right side.
// The trimmed string is returned both by value and reference.
//  in - str : ptr to char string to be trimmed
//  in - w   : which end(s) to trim:   both ends (B or b, the default),
//           :  left side (L or l), or right side (R or r).
//  out- str : ptr to the trimmed string by reference         
//  ret- str : ptr to the trimmed string by value
//-----------------------------------------------------------------------------
#ifdef NGS_PC_ENV
      char* trim_c(char* str, const char w);
#else

#  ifdef __cplusplus
    extern "C" {
#  endif

      char* trim_c(char* str, const char w);

#  ifdef __cplusplus
   }
#  endif

#endif //PC_ENV 

#endif // end TRIM_C_H
