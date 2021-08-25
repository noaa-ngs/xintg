#ifndef GRID_HEADER_H
#define GRID_HEADER_H

// %P%
// ----- constants ---------------------------------------------------
#pragma ident "$Id: grid_header.h 35301 2010-06-11 13:17:13Z Srinivas.Reddy $	20$Date: 2009/05/15 15:16:35 $ NGS"

// ----- standard library --------------------------------------------
// ----- classes, structs, types -------------------------------------
typedef struct {
    double lat_min;      // North Lat, decdeg (decimal degrees)
    double lon_min;      // East  Lon, decdeg
    double lat_delta;    //            decdeg
    double lon_delta;    //            decdeg
    long   lat_num;
    long   lon_num;
    long   ikind;
} GRID_HEADER;

#endif //~GRID_HEADER_H

