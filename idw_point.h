#ifndef IDW_POINT_H
#define IDW_POINT_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: idw_point.h 35301 2010-06-11 13:17:13Z Srinivas.Reddy $	20$Date: 2009/05/15 15:16:35 $ NGS"

// ----- standard library --------------------------------------------
// ----- classes, structs, types -------------------------------------
typedef struct {
    double latitude;       // North Lat, decdeg (decimal degrees)
    double longitude;      // East  Lon, decdeg
    double value;          // data value at position
} IDW_POINT;

#endif //~IDW_POINT_H

