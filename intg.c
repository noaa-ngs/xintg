// %P%
// ----- constants ---------------------------------------------------
static const char SCCSID[]="$Id: intg.c 118840 2020-10-19 15:16:45Z bruce.tran $	20$Date: 2010/06/21 16:47:51 $ NGS";
static const char PGMVER[]="3.18";
static const char PGMDAT[]="2012/12/13";
static const int  DEBUG = 0;           // diagnostics print if != 0
static const int  MEM_STEP = 40;       // dynamic allocation increment

// ----- standard library --------------------------------------------
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>      // sqrt()

// ----- classes, structures, types ----------------------------------
#include "grid_header.h"
#include "dataset1.h"

// ----- functions ---------------------------------------------------
#include "bb80ll.h"
#include "expform.h"
#include "c2v.h"
#include "ff1.h"
#include "ff1out.h"
#include "ff2.h"
#include "ff2out.h"
#include "ff4out.h"
#include "getdir_geoid.h"
//#include "which1.h"
//#include "getgrd_geoid.h"
#include "getgrd_vardis.h"
//#include "getheaders.h"
#include "intro.h"
//#include "interg.h"
#include "interg_idw.h"
#include "run_bbk.h"
#include "trim_c.h"


/*********************************************************************
* For further technical information, questions, or comments:
*   NOAA, National Geodetic Survey, N/NGS6, Silver Spring, MD  U.S.A.
*   Attn   : Daniel R. Roman, Ph.D.
*   Phone  : 301-713-3202
*   Fax    : 301-713-4172
*   e-mail : dan.roman@noaa.gov
*********************************************************************/

/*********************************************************************
*                                                                    *
*                            DISCLAIMER                              *
*                                                                    *
*   THIS PROGRAM AND SUPPORTING INFORMATION IS FURNISHED BY THE      *
* GOVERNMENT OF THE UNITED STATES OF AMERICA, AND IS ACCEPTED AND    *
* USED BY THE RECIPIENT WITH THE UNDERSTANDING THAT THE UNITED STATES*
* GOVERNMENT MAKES NO WARRANTIES, EXPRESS OR IMPLIED, CONCERNING THE *
* ACCURACY, COMPLETENESS, RELIABILITY, OR SUITABILITY OF THIS        *
* PROGRAM, OF ITS CONSTITUENT PARTS, OR OF ANY SUPPORTING DATA.      *
*                                                                    *
*   THE GOVERNMENT OF THE UNITED STATES OF AMERICA SHALL BE UNDER NO *
* LIABILITY WHATSOEVER RESULTING FROM ANY USE OF THIS PROGRAM.  THIS *
* PROGRAM SHOULD NOT BE RELIED UPON AS THE SOLE BASIS FOR SOLVING A  *
* PROBLEM WHOSE INCORRECT SOLUTION COULD RESULT IN INJURY TO PERSON  *
* OR PROPERTY.                                                       *
*                                                                    *
*   THIS PROGRAM IS PROPERTY OF THE GOVERNMENT OF THE UNITED STATES  *
* OF AMERICA.  THEREFORE, THE RECIPIENT FURTHER AGREES NOT TO ASSERT *
* PROPRIETARY RIGHTS THEREIN AND NOT TO REPRESENT THIS PROGRAM TO    *
* ANYONE AS BEING OTHER THAN A GOVERNMENT PROGRAM.                   *
*                                                                    *
*********************************************************************/

//Global variables
FILE *efp;
FILE *mfp;
FILE *nfp;
char old_cluster_rec[255];
int fatal_error;
int vc_unit;
int car97_unit;
int isNAD83 = 0;
int isPA11 = 0;
int isMA11 = 0;
int isFirstRun = 1;
char pcDir[255];
char hemisphere[2];
int imodel = 12; //starting model GEOID12B
#ifdef NGS_PC_ENV
int geoidCnt = 4;
#else
int geoidCnt = 16; //number of model to compute (GEOID12B, USGG2012, xGEOID14A, xGEOID14B, xGEOIDXXA, xGEOIDXXB)
#endif
int isFileFound = 1;
int is08_to_83 = 0; //use for testing consistency transformation from IGS08 to NAD83
int is_geoid12b = 1; //flag to indicate if processing geoid12b grids
double geoid_height12b = (double) -999.;
double ortho_height12b = 0.0;
int imodel12b = 12;
int count12b = 0; //use index 0 for geoid12b
char inputType[30]; //is input in NAD83 or IGS08
char outputType[30]; //IGS08 epoch 2022.00
char inputEpoch[13]; //2010.00 or 2005.00
char outputEpoch[13]; //2022.00
char outputVD[13]; //NAVD88
char newModel[13]; //GEOID18
char oldModel[13]; //GEOID12B
char currentModel[13];
//int isPR=0;
int isHTDPUsed = 0;
int isEast = 0;
int execHTDP = 0;
char saved_new_ellip_height[15];
char saved_new_latitude[50];
char saved_new_longitude[50];
double DN = 0.0;
double STD = 0.0;

void cat_date(char* prefix, char* cname);
int am_I_on_the_Pacific_plate(double fi, double la);

double interg(double xlat, double xlon, GRID_HEADER  vec_hdr[geoidCnt][50], FILE* vec_ifp[geoidCnt][50], int kk, int count);
double interg_idw(double xlat, double xlon,GRID_HEADER  vec_hdr[50], FILE* vec_ifp[50], int kk);
void getgrd_geoid(int imodel, char* dirnam, int is_subr, int* nfiles, int* nff, char vec_fnames[50][256], FILE* vec_ifp[geoidCnt][50], int count);

//local function
static double dms_dd(char *dms_str);
char * strtok_single (char * str, char const * delims);

int main( const int argc, const char* argv[] ) {
/*******************************************************************************
* NGS Program INTG - "INTerpolate Geoid"
* Interpolates the geoid height for a user specified position and geoid model.
* Estimates within gridded data models use spline or bilinear interpolation.
*
* Gridded data model files (*.bin) are direct access, unformatted, binary format.
* The order of bytes in the geoid model data files are
* --- depends on which platform the file was created ---
* Platform dependant endian condition is corrected for the binary data.
*******************************************************************************/
    FILE *fp;
    FILE *csv;
    int i = 0;
    double new_ellip_height[geoidCnt];
    double geoid_height[geoidCnt];
    double ortho_height[geoidCnt];
    char old_ellip_height[15];
    char station[50];
    char new_latitude[geoidCnt][50];
    char old_latitude[50];
    char new_longitude[geoidCnt][50];
    char old_longitude[50];
    char line[128];
    char ifname[256];
    char coords[6][50];
    char cinput[42];
    char cformat[42];
    char cformatFile[50];
    char modelA[15]; //model name selected
    char modelB[15];
    int heightA; //index array of each model height
    int heightB;
    int model;   //geoid model selected
    int latestModel = 20;
    int format; //selected output format
    char inputFrame[13];
    char outputFrame[13];
    char epochUsed[50];
    int moreLines = 0;
    double t = 2020.0;
    double t0 = 2020.0;

    for (i = 0; i <= geoidCnt; ++i) {
        strncpy (coords[i], "\0", sizeof(coords[i]));
    }

    for (i = 0; i < geoidCnt; ++i) {
        strncpy (new_latitude[i], "\0", sizeof(new_latitude[i]));
        strncpy (new_longitude[i], "\0", sizeof(new_longitude[i]));
    }
    strncpy(inputType, "\0", sizeof(inputType));
    strncpy(inputEpoch, "\0", sizeof(inputEpoch));
    strncpy(outputType, "\0", sizeof(outputType));
    strncpy(outputEpoch, "\0", sizeof(outputEpoch));
    strncpy(outputVD, "\0", sizeof(outputVD));
    strncpy(newModel, "\0", sizeof(newModel));
    strncpy(oldModel, "\0", sizeof(oldModel));
    strncpy(currentModel, "\0", sizeof(currentModel));
    strncpy(inputFrame, "\0", sizeof(inputFrame));
    strncpy(outputFrame, "\0", sizeof(outputFrame));
    strncpy(epochUsed, "\0", sizeof(epochUsed));
    strncpy (station, "\0", sizeof (station));
    strncpy (old_ellip_height, "\0", sizeof (old_ellip_height));
    strncpy (old_latitude, "\0", sizeof(old_latitude));
    strncpy (old_longitude, "\0", sizeof(old_longitude));
    strcpy(inputType,"IGS14 epoch 2010.00");
    strcpy(inputEpoch,"2010.00");
    strcpy(outputType,"IGS14 epoch 2020.00");
    strcpy(outputEpoch,"2020.00");
    strcpy(outputVD,"NAVD88");
    strcpy(newModel,"GEOID18");
    strcpy(oldModel,"GEOID12B");
    strcpy(inputFrame,"IGSI4");
    strcpy(outputFrame,"IGS14");
    strcpy(epochUsed,"IGS14(2010)  (EPOCH = 01-01-2020 (2020.00)");
    strncpy(cinput, "\0",  sizeof(cinput));
    strncpy(modelA, "\0",  sizeof(modelA));
    strncpy(modelB, "\0",  sizeof(modelB));

    //read in input file of coordinates and or ellipsoid height
    int iii = 0;
    while (iii == 0) {
#ifdef NGS_PC_ENV
        model = latestModel;
#else
        printf("\n");
        printf("\
            Which geoid model do you wish to use?\n\n\
               14 = xGEOID14\n\
               15 = xGEOID15\n\
               16 = xGEOID16\n\
               17 = xGEOID17\n\
               18 = xGEOID18\n\
               19 = xGEOID19\n\
               20 = xGEOID20\n\n");
        printf("Enter xGeoid model : ");
        strncpy(cinput, "\0", 42);
        fgets( cinput, 40, stdin);
        model = atoi(cinput);
#endif

        if ((model >= 14 && model <= latestModel)) {
            snprintf(modelA, sizeof (modelA), "xGEOID%d%s", model, "A");
            snprintf(modelB, sizeof (modelB), "xGEOID%d%s", model, "B");

            switch (model) {
            case 14:
                heightA = 2;
                heightB = 3;
                break;
            case 15:
                heightA = 4;
                heightB = 5;
                break;
            case 16:
                heightA = 6;
                heightB = 7;
                break;
            case 17:
                heightA = 8;
                heightB = 9;
                break;
            case 18:
                heightA = 10;
                heightB = 11;
                break;
            case 19:
                heightA = 12;
                heightB = 13;
                break;
            case 20:
                heightA = 14;
                heightB = 15;
                break;

            default:
                fprintf(stderr, "Error: Invalid option %d\n", model);
                exit(-1);
            }//~switch
#ifdef NGS_PC_ENV
                heightA = 2;
                heightB = 3;
#endif
            printf("\n");
        } else {
            fprintf(stderr,"Error: Not a valid response. Try again.\n");
            exit(-1);
        }
        printf("Enter your input file name : ");
        strncpy(ifname, "\0", sizeof(ifname));
        fgets(  ifname, 256, stdin);
        trim_c( ifname, 'b' );

        if ((fp = fopen( ifname, "r" )) == NULL) {
            printf("Error: Cannot find input file %s\nTry again\n", ifname);
            exit(1);
        } else {
            ++iii;
        }
        printf("\n");
        printf("\
            Please choose output format.\n\n\
               0 = JSON\n\
               1 = Screen\n\
               2 = File (csv format)\n\n");
        printf("Enter format : ");
        strncpy(cformat, "\0", 42);
        fgets( cformat, 40, stdin);
        format = atoi(cformat);
        if ((format == 0 || format == 1 || format == 2)) {
#ifdef NGS_PC_ENV
#else
            if (format == 2){
                printf("\n");
                printf("Enter output file: ");
                strncpy(cformatFile, "\0", 50);
                fgets( cformatFile, 48, stdin);
                trim_c( cformatFile, 'b' );
            }
#endif
        } else {
            fprintf(stderr, "Error: Invalid option %d\n", format);
            exit(-1);
        }
        printf("\n");
    }//~while


    strncpy(line, "\0", sizeof(line));
    int isFirstLine = 1;
    int len;
    int iRecord = 0;
    int isHeader = 1;
    int isCSV = 0;
    while ((int) fgets(line, 127, fp) != 0) {
        //printf ("out - %s", line);
        moreLines = 1;
        execHTDP = 0;
        strncpy (saved_new_latitude, "\0", sizeof(saved_new_latitude));
        strncpy (saved_new_longitude, "\0", sizeof(saved_new_longitude));

        // remove newline
        len = strlen(line);
        if( line[len-1] == '\n' ){
            line[len-1] = 0;
        }

        // remove ^M
        if( line[len-2] == 0x0d ){
            line[len-2] = 0;
        }

        if (isFirstLine){
            isFirstLine = 0;
            strncpy (hemisphere, "\0", sizeof(hemisphere));
            if(strstr(line, "NAD83") != NULL) {
                isNAD83 = 1;
                strcpy(inputType,"NAD83 (2011) epoch 2010.00");
                strcpy(inputEpoch,"2010.00");
                strcpy(inputFrame, "NAD83");
                strcpy(epochUsed,"NAD83(2011)  (EPOCH = 01-01-2010 (2010.00)");
            }
            if(strstr(line, "PA11") != NULL) {
                isPA11 = 1;
                strcpy(inputType,"NAD83 (PA11) epoch 2010.00");
                strcpy(inputEpoch,"2010.00");
                strcpy(inputFrame, "NAD83");
                strcpy(epochUsed,"NAD83(2011)  (EPOCH = 01-01-2010 (2010.00)");
            }
            if(strstr(line, "MA11") != NULL) {
                isMA11 = 1;
                strcpy(inputType,"NAD83 (MA11) epoch 2010.00");
                strcpy(inputEpoch,"2010.00");
                strcpy(inputFrame, "NAD83");
                strcpy(epochUsed,"NAD83(2011)  (EPOCH = 01-01-2010 (2010.00)");
            }


            if(strstr(line, "W") != NULL) {
                strcpy (hemisphere, "W");
            }else{
                strcpy (hemisphere, "E");
            }
//            printf ("\nYour input in %s\n",inputType);
//            printf("Your Result in IGS08: \n");
//            printf("----------------------------------------------------------------------------------\n");
        }else{
            char *token;
            int i = 0;

            token = strtok_single (line, ",");
            while (token) {
              //printf ("token - [%s]\n", *token ? token : "<empty>");
              strcpy(coords[i++], token);
              token = strtok_single (NULL, ",");
            }

            //printf("coords[0] [%s]\n",coords[0]);
            if (coords[0][0] == '\0'){
                exit(-1);
            }

            //printf("coors[3]: [%s]\n",coords[3]);
            if (strcmp(coords[3], "") != NULL){
                t = atof(coords[3]);
            }else{
                t = 2020.0f;
            }

            //printf("model[%d]\n",model);
            if (model == 19){
                is_geoid12b = 0; //process non-geoid12b grids
                process (24, 2, &coords, &new_latitude, &new_longitude, &new_ellip_height, &geoid_height, &ortho_height);
                process (25, 3, &coords, &new_latitude, &new_longitude, &new_ellip_height, &geoid_height, &ortho_height);
            }else{
            //compute GEOID12B separately
            //printf("compute GEOID18\n");
            strcpy(currentModel, newModel);
            is_geoid12b = 1; //process geoid12b grids
            process (imodel12b, count12b, &coords, &new_latitude, &new_longitude, &new_ellip_height, &geoid_height, &ortho_height);

            //if no geoid height exists for GEOID18
            //then compute GEOID12B
            //printf("geoid_height[0][%lf]\n",geoid_height[0]);
            if (geoid_height[0] == -999.) {
                strcpy(currentModel, oldModel);
                process (100, count12b, &coords, &new_latitude, &new_longitude, &new_ellip_height, &geoid_height, &ortho_height);
            }
            //printf("done GEOID18\n");
#ifdef NGS_PC_ENV
            //printf("compute 20A/20B\n");
            is_geoid12b = 0; //process non-geoid12b grids
            //process (25, 1, &coords, &new_latitude, &new_longitude, &new_ellip_height, &geoid_height, &ortho_height);
            process (26, 2, &coords, &new_latitude, &new_longitude, &new_ellip_height, &geoid_height, &ortho_height);
            //printf("geoid_height20A[%lf]\n",geoid_height[2]);
            //printf("compute 20B\n");
            process (27, 3, &coords, &new_latitude, &new_longitude, &new_ellip_height, &geoid_height, &ortho_height);
            //printf("geoid_height20B[%lf]\n",geoid_height[3]);
#else
            //printf("increment imodel[%d]\n",imodel);
            imodel = imodel + 1; //increment 1 since geoid12b already processed

            //compute xGEOID separately
            //printf("geoidCnt[%d]\n",geoidCnt);
            for (i = 1; i < geoidCnt; ++i) {
                is_geoid12b = 0; //process non-geoid12b grids
                //printf("process count[%d]\n",i);
                process (imodel++, i, &coords, &new_latitude, &new_longitude, &new_ellip_height, &geoid_height, &ortho_height);
            }
#endif
            }


            //if user entered DMS return DMS
            if (strstr(coords[0], " ") != NULL){
                //do nothing
            }else{//if user entered decimal degree return decimal degree
                char tc[50];
                double c;

                c = dms_dd (new_latitude[2]);
                sprintf(tc,"%13.10f",c);
                strcpy (new_latitude[1], tc);

                c = dms_dd (new_longitude[2]);
                sprintf(tc,"%13.10f",c);
                strcpy (new_longitude[1], tc);

                new_ellip_height[1] = new_ellip_height[2];

            }

            if (isEast && isHTDPUsed){
                double n = atof(new_longitude[1]);
                n = 360. - n;
                char tc[50];
                sprintf(tc,"%13.10f",n);
                strcpy (new_longitude[1], tc);
            }

            //printf("model[%d]\n",model);
            if (model == 19){
                if (geoid_height[2] == -999. ) {
                    printf("Point coordinate outside the area for which the model is valid\n");
                    printf("If area is Alaska, Hawaii, Guam, American Samoa, or Commonwealth of Northern Marianas Islands please use GEOID12B\n");
                }else {
                    if (format == 0) {
                        if (isHeader) {
                            printf("{");
                            printf("%s:%s,", "\"n\"", "\"N: the geoid height at epoch t0 = 2020.0, which is geocentric and relative to the GRS80 reference ellipsoid.\"");
                            printf("%s:%s,", "\"dn\"", "\"DN: the time-dependent geoid change computed between user inputted epoch (t) and t0.\"");
                            printf("%s:%s,", "\"add\"", "\"To obtain the geoid height at user inputted epoch (t), add N + DN.  Either Model A or Model B N values may be used for this depending on user preference.\"");
                            printf("%s:%s,", "\"link\"", "\"For the orthometric height change from NAVD88 to NAP2022, please visit https://www.ngs.noaa.gov/datums/newdatums/WhatToExpect.shtml.\"");
                            printf("%s:%s", "\"inputs\"", "[");
                            printf("{%s:%d,", "\"Cnt\"", iRecord++);
                            isHeader = 0;
                        } else {
                            printf(",{%s:%d,", "\"Cnt\"", iRecord++);
                        }

                        if (moreLines) {
                            moreLines = 0;

                            printf("%s:\"%s\",", "\"Station\"", coords[4]);
                            printf("\"Latitude\":\"%s\",", coords[0]);
                            printf("\"Longitude\":\"%s\",", coords[1]);
                            printf("\"%s_Ht\":\"%.3lf\",", modelA, geoid_height[2]);
                            printf("\"%s_Ht\":\"%.3lf\",", modelB, geoid_height[3]);
                            printf("\"DN\":\"%.3lf\",", (t - t0) * DN);
                            printf("\"Epoch\":\"%.1f\"",t);
                            printf("}");

                        }

                    }else if (format == 1) {
                        if (isHeader) {
                            printf("N: the geoid height at epoch t0 = 2020.0, which is geocentric and relative to the GRS80 reference ellipsoid.\n");
                            printf("DN: the time-dependent geoid change computed between user inputted epoch (t) and t0.\n");
                            printf("To obtain the geoid height at user inputted epoch (t), add N + DN.  Either Model A or Model B N values may be used for this depending on user preference.\n");
                            printf("For the orthometric height change from NAVD88 to NAP2022, please visit https://www.ngs.noaa.gov/datums/newdatums/WhatToExpect.shtml.\n");
                            printf("\n");
                            printf("%-30s Latitude        Longitude        N (Model A)  N (Model B)  DN       Epoch\n", " ");
                            printf("Station Name %-25s                          meters       meters       meters\n", " ");
                            isHeader = 0;
                        }

                        //printf("t: %lf\n",t);
                        //printf("t0: %lf\n",t0);
                        //printf("DN: %lf\n",DN);
                        printf("%-30s %-15s %-15s %-8.3lf     %-8.3lf   %8.3lf    %-6.1lf\n", coords[4], coords[0], coords[1], geoid_height[2], geoid_height[3], (t - t0) * DN, t);
                        //printf("\n");
                    }else {
                        //writing to file
                        //header columns
                        if (isHeader) {
                            //create temp output file
                            char cfile[256];
                            char datetime[256];
                            char csv_out[256];
                            char dummy[2];

                            strncpy(datetime, "\0", sizeof (datetime));
                            strncpy(dummy, "\0", sizeof (dummy));
                            strncpy(cfile, "\0", sizeof (cfile));
                            strncpy(csv_out, "\0", sizeof (csv_out));


#ifdef NGS_PC_ENV
                            strcat(cfile, pcDir);
                            strcpy(cfile, "geoidout.csv");

                            //check if file exist and delete
                            if (access(cfile, F_OK) != -1) {
                                // file exists
                                remove(cfile);
                            } else {
                                // file doesn't exist
                            }
#else
                            strcat(cfile, cformatFile);
#endif

                            if ((csv = fopen(cfile, "a")) == NULL) {
                                printf("\n ABORT: Can not open %s.\n", cfile);
#ifndef NGS_PC_ENV
                                exit(-1);
#else
                                exit(1);
#endif
                            }

                            isHeader = 0;
                            isCSV = 1;
                            fprintf(csv, "N: the geoid height at epoch t0 = 2020.0, which is geocentric and relative to the GRS80 reference ellipsoid.\n");
                            fprintf(csv, "DN: the time-dependent geoid change computed between user inputted epoch (t) and t0.\n");
                            fprintf(csv, "To obtain the geoid height at user inputted epoch (t), add N + DN.  Either Model A or Model B N values may be used for this depending on user preference.\n");
                            fprintf(csv, "For the orthometric height change from NAVD88 to NAP2022, please visit https://www.ngs.noaa.gov/datums/newdatums/WhatToExpect.shtml.\n");
                            fprintf(csv, "Cnt,Station,Latitude,Longitude,%s_Ht,%s_Ht,DN,Epoch\n",modelA, modelB);
                        }

                        //values of each input record
                        fprintf(csv, "%d,%s,%s,%s,%.3lf,%.3lf,%.3lf,%.1lf\n",
                                iRecord++, coords[4], coords[0], coords[1], geoid_height[2], geoid_height[3], (t - t0) * DN, t);

                        fflush(csv);

                    }


                }
            }else{
                //0 - JSON; 1 - screen ; 2 - csv file
                if (format == 0){
                    if (isHeader){
                        printf("{");
                        printf("%s:%s,", "\"note\"", "\"Note: The GRS80 ellipsoid is used for both NAD83 and IGS14.\"");
                        printf("%s:%s,", "\"n\"", "\"N: the geoid height at epoch t0 = 2020.0, which is geocentric and relative to the GRS80 reference ellipsoid.\"");
                        printf("%s:%s,", "\"accuracy\"", "\"Accuracy: Estimated 95%% confidence interval for geoid height.\"");
                        printf("%s:%s,", "\"dn\"", "\"DN: the time-dependent geoid change computed between user inputted epoch (t) and t0. To obtain the dynamic geoid height at user inputted epoch (t), add N + DN.  Either Model A or Model B N values may be used for this depending on user preference.\"");
                        printf("%s:%s","\"inputs\"","[");
                        printf("{%s:%d,","\"Cnt\"",iRecord++);
                        isHeader = 0;
                    }else{
                        printf(",{%s:%d,","\"Cnt\"",iRecord++);
                    }

                    if (moreLines){
                        moreLines = 0;

                        printf("%s:\"%s\",","\"Station\"",coords[4]);
                        printf("\"%s_Lat\":\"%s\",",inputFrame,coords[0]);
                        printf("\"%s_Lon\":\"%s\",",inputFrame,coords[1]);
                        printf("\"%s_Eht\":\"%s\",",inputFrame,coords[2]);
                        printf("%s:\"%s\",","\"Input_Epoch\"",inputEpoch);
                        printf("\"%s_Lat\":\"%s\",",outputFrame,new_latitude[1]);
                        printf("\"%s_Lon\":\"%s\",",outputFrame,new_longitude[1]);
                        printf("\"%s_Eht\":%.3lf,",outputFrame,new_ellip_height[1]);
                        printf("%s:\"%s\",","\"Output_Epoch\"",outputEpoch);
                        printf("%s:%.3lf,","\"GEOID18_Ht\"",geoid_height[0]);
                        printf("%s%s%s:%.3lf,","\"Oht_",outputVD,"\"",ortho_height[0]);
                        printf("\"%s_Ht\":%.3lf,",modelA,geoid_height[heightA]);
                        printf("\"%s_Ht\":%.3lf,",modelB,geoid_height[heightB]);
                        printf("\"%s_Accuracy\":%.3lf,",modelB,STD);
                        //printf("%s:%.3lf,","\"Oht_USGG2012\"",ortho_height[1]);
                        printf("\"Oht_%s\":%.3lf,",modelA,ortho_height[heightA]);
                        printf("\"Oht_%s\":%.3lf,",modelB,ortho_height[heightB]);
                        //printf("%s%s%s:%.3lf,","\"Oht_Diff(USGG2012-",outputVD,")\"",ortho_height[1] - ortho_height[0]);
                        printf("\"Oht_Diff(%s-%s)\":%.3lf,",modelA,outputVD,ortho_height[heightA] - ortho_height[0]);
                        printf("\"Oht_Diff(%s-%s)\":%.3lf,",modelB,outputVD,ortho_height[heightB] - ortho_height[0]);
                        printf("\"DN\":%.3lf,",(t-t0) * DN);
                        printf("\"Epoch\":%.3lf",t);
                        printf("}");

                    }

                }else if (format == 1){
					//printf("format == 1\n");
                    if (isHeader){
                        printf("Note: The GRS80 ellipsoid is used for both NAD83 and IGS14.\n");
                        printf("N: the geoid height at epoch t0 = 2020.0, which is geocentric and relative to the GRS80 reference ellipsoid.\n");
                        printf("Accuracy: Estimated 95%% confidence interval for geoid height.\n");
                        printf("DN: the time-dependent geoid change computed between user inputted epoch (t) and t0. To obtain the dynamic geoid height at user inputted epoch (t), add N + DN.  Either Model A or Model B N values may be used for this depending on user preference.\n\n");
                        isHeader = 0;
                    }
                    //printf("geoid_height[0][%lf]\n",geoid_height[0]);
                    if (geoid_height[0] != -999.) {
                        if (strcmp(coords[2], " ") == NULL) {
                            printf("\nYour input in %s:\n", inputType);
                            printf("Latitude %8sLongitude %8sEllipsoid Height%8sStation\n", " ", " ", " ");
                            printf("%-15s  %-15s   %-20s    %s\n", coords[0], coords[1], " ", coords[4]);
                            printf("\nYour Result in %s: \n", outputType);
                            printf("Latitude %8sLongitude %8sEllipsoid Height\n", " ", " ");
                            printf("%-15s  %-15s   %-20s\n", new_latitude[1], new_longitude[1], " ");
                            printf("\n");
                            printf("Geoid Model %8sGeoid Height(m) %8s Ortho(model)-Ortho(%s)(m)\n", " ", " ", currentModel);
                            //printf ("GEOID12B  %24.3lf %25.3lf\n", geoid_height[0], geoid_height[heightB]-geoid_height[0]);
                            //printf("xGEOID19B  %24.3lf %25.3lf\n", geoid_height[1], geoid_height[heightB] - geoid_height[1]);
                            printf("%s %24.3lf %25.3lf\n", modelA, geoid_height[heightA], geoid_height[heightB] - geoid_height[heightA]);
                            printf("%s %24.3lf %25.3lf\n", modelB, geoid_height[heightB], geoid_height[heightB] - geoid_height[heightB]);
                            //printf("*Changes in the orthometric height relative to NAVD88(GEOID12B)\n");
                            printf("--------------------------------------------------------------------------------------------\n");
                        } else {
                            printf("\nYour input in %s:\n", inputType);
                            printf("Latitude %8sLongitude %8sEllipsoid Height%8sStation\n", " ", " ", " ");
                            printf("%-15s  %-15s   %-20s    %s\n", coords[0], coords[1], coords[2], coords[4]);
                            printf("\nYour Result in %s: \n", outputType);
                            printf("Latitude %8sLongitude %8sEllipsoid Height\n", " ", " ");
                            if (strcmp(coords[2], "") == NULL) {
								printf("%-15s  %-15s   %-20s\n", new_latitude[1], new_longitude[1], " ");
							}else{
								printf("%-15s  %-15s   %-20.3lf\n", new_latitude[1], new_longitude[1], new_ellip_height[1]);
							}

                            printf("\n");
                            printf("The geoid height of %s (with respect to NAD83): %17.3lf m\n", currentModel, geoid_height[0]);
                            if (strcmp(coords[2], "") == NULL) {
                                printf("The orthometric height in %s (based on %s): %17s     \n", outputVD, currentModel, " ");
                            }else{
                                printf("The orthometric height in %s (based on %s): %17.3lf m\n", outputVD, currentModel, ortho_height[0]);
                            }
                            printf("\n");
                            printf("Estimated orthometric height in North American-Pacific Geopotential Datum of 2022 (NAPGD2022) \n");
                            printf("based on different geoid models (all quantities in meters):\n");
                            printf("\n");
                            printf("Geoid Model %3sGeoid Height %3sAccuracy %3sOrtho Height %3sOrtho(model)-%s(%s) %3sDN%5sEpoch\n", " ", " ", " ", " ", outputVD, currentModel, " ", " ");
                            //printf ("GEOID12B  %24.3lf %24.3lf %25.3lf\n", geoid_height[0], ortho_height[0], ortho_height[heightB]- ortho_height[0]);
                            //printf("xGEOID19B %14.3lf %13s %14.3lf %21.3lf %s\n", geoid_height[1], " ", ortho_height[1], ortho_height[1] - ortho_height[0], " ");
                            if (strcmp(coords[2], "") == NULL) {
                                printf("%s %14.3lf %13s %14s %21s %19.3lf %8.1lf\n", modelA, geoid_height[heightA], " ", " ", " ", (t - t0) * DN, t);
                                printf("%s %14.3lf %12.3lf %15s %21s %19.3lf %8.1lf\n", modelB, geoid_height[heightB], STD, " ", " ", (t - t0) * DN, t);
                            }else{
                                printf("%s %14.3lf %13s %14.3lf %21.3lf %19.3lf %8.1lf\n", modelA, geoid_height[heightA], " ", ortho_height[heightA], ortho_height[heightA] - ortho_height[0], (t - t0) * DN, t);
                                printf("%s %14.3lf %12.3lf %15.3lf %21.3lf %19.3lf %8.1lf\n", modelB, geoid_height[heightB], STD, ortho_height[heightB], ortho_height[heightB] - ortho_height[0], (t - t0) * DN, t);
                            }
                            //printf("*Changes in the orthometric height relative to NAVD88(GEOID12B)\n");
                            printf("--------------------------------------------------------------------------------------------\n");
                        }
                    } else {
                        if (strcmp(coords[2], " ") == NULL) {
                            printf("\nYour input in %s:\n", inputType);
                            printf("Latitude %8sLongitude %8sEllipsoid Height%8sStation\n", " ", " ", " ");
                            printf("%-15s  %-15s   %-20s    %s\n", coords[0], coords[1], " ", coords[4]);
                            printf("\nYour Result in %s: \n", outputType);
                            printf("Latitude %8sLongitude %8sEllipsoid Height\n", " ", " ");
                            printf("%-15s  %-15s   %-20s\n", new_latitude[1], new_longitude[1], " ");
                            printf("\n");
                        } else {
                            printf("\nYour input in %s:\n", inputType);
                            printf("Latitude %8sLongitude %8sEllipsoid Height%8sStation\n", " ", " ", " ");
                            printf("%-15s  %-15s   %-20s    %s\n", coords[0], coords[1], coords[2], coords[4]);
                            printf("\nYour Result in %s: \n", outputType);
                            printf("Latitude %8sLongitude %8sEllipsoid Height\n", " ", " ");
                            printf("%-15s  %-15s   %-20.3lf\n", new_latitude[1], new_longitude[1], new_ellip_height[1]);
                            printf("\n");
                        }
                        printf("One or more errors have occurred.\n");
                        printf("\n");
                        printf("Errors may be due to invalid input, such as\n");
                        printf("\n");
                        printf("Improperly formatted input latitude or longitude\n");
                        printf("Point coordinate outside the area for which the model is valid\n");
                        printf("\n");
                        printf("Point is out of bounds.\n");
                        printf("--------------------------------------------------------------------------------------------\n");
                    }
                }else{
                    //writing to file
                    //header columns
                    if (isHeader){
                        //create temp output file
                        char cfile[256];
                        char datetime[256];
                        char csv_out[256];
                        char dummy[2];

                        strncpy(datetime, "\0", sizeof (datetime));
                        strncpy(dummy, "\0", sizeof (dummy));
                        strncpy(cfile, "\0", sizeof (cfile));
                        strncpy(csv_out, "\0", sizeof (csv_out));


    #ifdef NGS_PC_ENV
                        strcat(cfile, pcDir);
                        strcpy(cfile, "geoidout.csv");

                        //check if file exist and delete
                        if( access( cfile, F_OK ) != -1 ) {
                                                // file exists
                                                remove(cfile);
                                            } else {
                                                // file doesn't exist
                                            }
    #else
                        strcat (cfile, cformatFile);
    #endif

                        if ((csv = fopen(cfile, "a")) == NULL) {
                            printf("\n ABORT: Can not open %s.\n", cfile);
    #ifndef NGS_PC_ENV
                            exit(-1);
    #else
                            exit(1);
    #endif
                        }

                        isHeader = 0;
                        isCSV = 1;
                        fprintf(csv,"Note: The GRS80 ellipsoid is used for both NAD83 and IGS14.\n");
                        fprintf(csv,"N: the geoid height at epoch t0 = 2020.0, which is geocentric and relative to the GRS80 reference ellipsoid.\n");
                        fprintf(csv,"Accuracy: Estimated 95%% confidence interval for geoid height.\n");
                        fprintf(csv,"DN: the time-dependent geoid change computed between user inputted epoch (t) and t0. To obtain the dynamic geoid height at user inputted epoch (t), add N + DN.  Either Model A or Model B N values may be used for this depending on user preference.\n\n");

                        fprintf(csv,"Cnt,Station,%s_Lat,%s_Lon,%s_Eht,Input_Epoch,%s_Lat,%s_Lon,%s_Eht,Output_Epoch,GEOID18_Ht,Oht_%s,"
                                "%s_Ht,%s_Ht,%s_Accuracy,Oht_%s,Oht_%s,"
                                "Oht_Diff(%s-%s),Oht_Diff(%s-%s),DN,Epoch\n",
                                inputFrame,inputFrame,inputFrame,outputFrame,outputFrame,outputFrame,outputVD,
                                modelA,modelB,modelB,modelA,modelB,outputVD,modelA,outputVD,modelB,outputVD);
                    }

                    //values of each input record
                    fprintf(csv,"%d,%s,%s,%s,%s,%s,%s,%s,%.3lf,%s,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf\n",
                            iRecord++,coords[4],coords[0],coords[1],coords[2],inputEpoch,
                            new_latitude[1],new_longitude[1],new_ellip_height[1],outputEpoch,geoid_height[0],ortho_height[0],
                            geoid_height[heightA],geoid_height[heightB],STD,ortho_height[heightA],ortho_height[heightB],
                            ortho_height[heightA] - ortho_height[0],ortho_height[heightB] - ortho_height[0],(t - t0) * DN,t);

                    fflush(csv);

                }
            }


            //clear coords
            for (i = 0; i <= geoidCnt; ++i) {
                strncpy (coords[i], "\0", sizeof(coords[i]));
            }

            //reset imodel
            imodel = 12;
        }
    }

    if (format == 0 && moreLines == 0){
        printf("]}\n");
    }

    //if csv file was open then close it
    if (isCSV == 1){
#ifdef NGS_PC_ENV
		printf ("\ncsv output file - %s%s\n",pcDir,"geoidout.csv");
#endif
        fclose(csv);
    }
    fclose(fp);

}//~intg

double dms_dd(char *dms_str)
{

   int deg = 0;
   int min = 0;

   double sec  = 0.0;
   double dd   = 0.0;

   char* err_ptr;

   char deg_str[4];
   char min_str[3];
   char sec_str[12];
   char coord[3][50];
   char fn[] = "dms_dd()";

   strncpy(deg_str, "\0",  4);
   strncpy(min_str, "\0",  4);
   strncpy(sec_str, "\0", 12);

   if( strlen(dms_str) > 17 ) {
      printf("String into %s too long\n", fn);
      printf("Max length := dddmmss.12345678\n");
      printf("String      = %s\n", dms_str );
      exit(1);
   }

     char *token;
     int i = 0;
     token = strtok(dms_str, " ");
     do
     {
       //printf("token: %s\n", token);
       strcpy(coord[i++], token);
     }
     while (token = strtok(NULL, " "));

    strcpy(deg_str, coord[0]);
    strcpy(min_str, coord[1]);
    strcpy(sec_str, coord[2]);


   deg = strtol(deg_str, &err_ptr, 10);
   //printf ("%f\n",deg);
   if (strlen(err_ptr) > 0) {
      printf("ERROR: Invalid degrees in %s: %s\n", fn, err_ptr );
      printf("       DMS entered = %s\n", dms_str );
      exit(1);
   }

   min = strtol(min_str, &err_ptr, 10);
   //printf ("%f\n",min);
   if (strlen(err_ptr) > 0) {
      printf("ERROR: Invalid minutes in %s: %s\n", fn, err_ptr );
      printf("       DMS entered = %s\n", dms_str );
      exit(1);
   }

   sec = strtod(sec_str, &err_ptr);
   if (strlen(err_ptr) > 0) {
      printf("ERROR: Invalid seconds in %s: %s\n", fn, err_ptr );
      printf("       DMS entered = %s\n", dms_str );
      exit(1);
   }

   //check if negative degree then make sure min/sec are negative
   if (deg < 0){
       dd = deg - min/60.0 - sec/3600.0;    // decimal degrees
   }else{
       dd = deg + min/60.0 + sec/3600.0;    // decimal degrees
   }

   return( dd );

}//~dms_dd



/**
 * process inputs from user
 * @param imodel model number being called, in order USGG2012, GEOID12B, xGEOID14A, xGEOID14B, xGEOID15A, xGEOID15B
 * @param count number of calls used for array
 * @param coords coordinates array of lat, lon, and eht
 * @param new_latitude new latitude if converted
 * @param new_longitude new longitude if converted
 * @param new_ellip_height new ellipsoid height if converted
 * @param geoid_height geoid heights array returned
 * @param ortho_height orthometric heights array returned, if any
 * @return
 */
int process (int imodel, int count, char coords[geoidCnt][50], char new_latitude[geoidCnt][50], char new_longitude[geoidCnt][50], double new_ellip_height[], double geoid_height[], double ortho_height[]){
    FILE* ifp;
    FILE* ofp;
    FILE* vec_ifp[geoidCnt][50];            // vector of FILE* of height grids
    FILE* var_ifp[50];            // vector of FILE* of variance grids
    FILE* vdar_ifp[50];            // vector of FILE* of variance grids
    //FILE* dis_ifp[50];            // vector of FILE* of distance grids
    char  vec_fnames[50][256];    // vector of filenames of height grids
    char  var_fnames[50][256];    // vector of filenames of variance grids
    char  vdar_fnames[50][256];    // vector of filenames of variance grids
    //char  dis_fnames[50][256];    // vector of filenames of distance grids

    GRID_HEADER vec_hdr[geoidCnt][50];      // vector of header file data
    GRID_HEADER var_hdr[50];      // vector of header file data
    GRID_HEADER vdar_hdr[50];      // vector of header file data
    //GRID_HEADER dis_hdr[50];      // vector of header file data
    DATASET1*   vec_data;  // ptr to vector of lat-lon-text input data

    char    dash70[]="----------------------------------------------------------------------";
    char    card[90];
    char    card2[90];
    char    card86[90];
    char    dirnam[256];
    char    ddirnam[256];
    char    rec_in[90];
    char    ifname[256];  //  input file name
    char    ofname[256];  // output file name
    char    ofyn[2];
    char    cinput[42];
    char    text[42];
    double  dh[geoidCnt];

    double  xlat;
    double  xlon;

    // Variables for statistics of the run
    double  minght =  999999.0; // smallest geoid value
    double  minlat =      90.0;    // lat at min location
    double  minlon =     360.0;    // lon at min location
    double  maxght = -999999.0; // largest  geoid value
    double  maxlat =     -90.0;    // lat at max location
    double  maxlon =    -180.0;    // lon at max location
    double  xn     = -999999.0;    // northernmost lat
    double  xs     =  999990.0;   // southernmost lat
    double  xe     = -999999.0;    // easternmost lon
    double  xw     =  999999.0;    // westernmost lon
    double  fact   = 0.0;
    double  ave    = 0.0;
    double  std    = 0.0;
    double  rms    = 0.0;
    double  geoidHt= 0.0;         // solution
    double  variance = -999.;
    double  vdariance = -999.;
    double  stddev = -999.;
    double  distance = -999.;

    int  iform  = 0;
    int  ipos;
    int  keep   = 0;
    int  kount  = 0;
    int  mem_limit = 0;
    int  kk  = 0;
    int  mm  = 0;
    int  iii = 0;
    int  iinput;
    int  is_subr;
    int  nfiles = 0, nvfiles = 0, nvdfiles = 0;
    int  nff, nvff, ndff;
    int  ii;
    int  jj;
    int  poseast = 0;  // Assume west longitude
    int  outfil  = 0;  // Cstd:  0 = false;  !0 = not false
    char old_ellip_height[15];
    char station[50];
    char old_latitude[50];
    char old_latitude2[50];
    char old_longitude[50];
    char old_longitude2[50];


    // -----------------------------------------------------
    // Initialize variables - Zero out the statistics for this run
    // -----------------------------------------------------
    strncpy(dirnam, "\0", sizeof (dirnam));
    strncpy(ddirnam, "\0", sizeof (ddirnam));
    strncpy(cinput, "\0",  sizeof (cinput));
    strncpy(text,   "\0",  sizeof (text));
    strncpy(card,   "\0",  sizeof (card));
    strncpy(card2,  "\0",  sizeof (card2));
    strncpy(card86, "\0",  sizeof (card86));
    strncpy(rec_in, "\0",  sizeof (rec_in));
    strncpy(ofyn,   "\0",   sizeof (ofyn));
    strncpy (station, "\0", sizeof (station));
    strncpy (old_ellip_height, "\0", sizeof (old_ellip_height));
    strncpy (old_latitude, "\0", sizeof(old_latitude));
    strncpy (old_latitude2, "\0", sizeof(old_latitude2));
    strncpy (old_longitude, "\0", sizeof(old_longitude));
    strncpy (old_longitude2, "\0", sizeof(old_longitude2));


    for (ii = 0; ii < 50; ++ii) {
        strncpy(vec_fnames[ii], "\0", sizeof(vec_fnames[ii]));
        strncpy(var_fnames[ii], "\0", sizeof(var_fnames[ii]));
        strncpy(vdar_fnames[ii], "\0", sizeof(vdar_fnames[ii]));
        //strncpy(dis_fnames[ii], "\0", 256);
    }


    iii = 0;

    //get values from coords
    if (hemisphere[0] == 'E'){
        poseast = 1;
    }

    strcpy (old_latitude, coords[0]);
    strcpy (old_latitude2, coords[0]); //temp
    strcpy (old_longitude, coords[1]);
    strcpy (old_longitude2, coords[1]); //temp
    strcpy (old_ellip_height, coords[2]);
    strcpy (station, coords[4]);

    // ---------------------------------------------------------
    // Which directory are the geoid files in?
    // For Unix platform, set the directory location for the user.
    // ---------------------------------------------------------
    getdir_geoid(imodel, dirnam);

    if (imodel >= 24 && imodel <= 27) {
        getdir_dgeoid(imodel, ddirnam);
    }

    // ---------------------------------------------------------
    // Create the list of files that must be opened, and open them.
    // Return with a count of how many files were opened,
    // and a flag (ios) of which files are open and which are not.
    // IS_SUBR : run as subroutine: false=0; true=(not zero) (c std notation)
    // ---------------------------------------------------------
    is_subr = 0;
    int dummyCount = 0;
    if ( (count < geoidCnt)){
        if (count == geoidCnt -1){
            isFileFound = 0;
        }
        getgrd_geoid(imodel, dirnam, is_subr, &nfiles, &nff, vec_fnames, vec_ifp, count);
    }


    if (imodel >= 24 && imodel <= 27) {
        getgrd_dgeoid(imodel, ddirnam, is_subr, &nvdfiles, &ndff, vdar_fnames, vdar_ifp, dummyCount);
        getheaders( vdar_ifp, vdar_hdr, nvdfiles, dummyCount );
    }
    // ---------------------------------------------------------
    // Read the headers of all geoid files which
    // where opened, and store that information.
    // Compute the max lat/lon from the header information.
    // Apply endian correction if required
    // ---------------------------------------------------------
    getheaders( vec_ifp, vec_hdr, nfiles, count );

    if (imodel == 27){
        getgrd_vardis(imodel, dirnam, is_subr, &nvfiles, &nvff, &ndff, var_fnames, var_ifp);
        getheaders( var_ifp, var_hdr, nvfiles, dummyCount );
    }

    //getheaders( dis_ifp, dis_hdr, nvdfiles );



    // ======================================================================
    // Now handle the 4 input types
    // ======================================================================

    mem_limit += MEM_STEP;
    vec_data = (DATASET1*) calloc(mem_limit, sizeof(DATASET1));
    if (vec_data == NULL ) {
        fprintf(stderr, "Out of system memory - allocation fails\n");
#ifndef NGS_PC_ENV
        exit(-1);
#endif
    }

    // -----------------------------------------------
    // Input by file, free format type 1
    // -----------------------------------------------
    if (iform == 1) {
        while( fgets(rec_in, 90, ifp) ) {
            trim_c(rec_in, 'r');
            if (strlen(rec_in) < 5) continue;

            // Find the lat/lon value
            // Longitude always returns as 0->360...whether this
            // is positive east or west is fixed a few lines down
            strncpy(text, "\0", 42);
            ff1(rec_in, &xlat, &xlon, text);

            // If the lat/lon values came back as -999, set the
            // geoid value to -999 and skip the interpolation
            if (xlat == -999. || xlon == -999.) {
                geoidHt = (double)-999.;
                continue;
            }

            // Force input West longitude to be positive East
            if (poseast == 0) {           // C std: 0 = false; is NOT east lon
                xlon = 360. - xlon;
            }
            // Now have the lat/lon pair from input file type 2

            DATASET1 thisSet;
            strncpy(thisSet.text, "\0", 42);
            thisSet.lat     = xlat;
            thisSet.lon     = xlon;
            thisSet.poseast = poseast;
            strcpy(thisSet.text, text);
            vec_data[kount] = thisSet;
            ++kount;
            if (kount >= mem_limit) {
                mem_limit += MEM_STEP;
            vec_data = (DATASET1*)realloc(vec_data, mem_limit*sizeof(DATASET1));
                if (vec_data == NULL) {
                   printf("Out of system memory - allocation fails\n");
                   exit(-1);
                }
            }
        }//~while
    }//~if (iform = 1)

    // -----------------------------------------------
    // Input file, free format type 2
    // -----------------------------------------------
    else if (iform == 2) {
        while( fgets(rec_in, 90, ifp) ) {
            trim_c(rec_in, 'r');

            // Find the lat/lon value
            // Longitude always returns as 0->360...whether this
            // is positive east or west is fixed a few lines down
            strncpy(text, "\0", 42);
            ff2(rec_in, &xlat, &xlon, text);

            // If the lat/lon values came back as -999, set the
            // geoid value to -999 and skip the interpolation
            if (xlat == -999. || xlon == -999.) {
                geoidHt = (double) -999.;
                continue;
            }

            // Force input West longitude to be positive East
            if (poseast == 0) {           // C std: 0 = false; is NOT east lon
                xlon = 360. - xlon;
            }
            // Now have the lat/lon pair from input file type 2

            DATASET1 thisSet;
            strncpy(thisSet.text, "\0", 50);
            thisSet.lat     = xlat;
            thisSet.lon     = xlon;
            thisSet.poseast = poseast;
            strcpy(thisSet.text, text);
            vec_data[kount] = thisSet;
            ++kount;
            if (kount >= mem_limit) {
                mem_limit += MEM_STEP;
            vec_data = (DATASET1*)realloc(vec_data, mem_limit*sizeof(DATASET1));
                if (vec_data == NULL) {
                    printf("Out of system memory - allocation fails\n");
                    exit(-1);
                }
            }

        }//~while
    }//~if (iform = 2)

    // -----------------------------------------------
    // Input by file, horizontal bluebook
    // Function contains all processes.
    // -----------------------------------------------
    else if (iform == 3) {
        //run_bbk( ifp, ofp, vec_ifp, vec_hdr, vec_fnames, nff, imodel );

        //fclose(ifp);
        //fclose(ofp);
        return(0);

    }//~(iform == 3) bluebook

    // -----------------------------------------------
    // Input by prompts
    // -----------------------------------------------
    else {
        //jj = 0;
        //do {
            xlat = c2v(old_latitude2, 1);
            if (fabs(xlat + 999) < 0.001) {
                printf("Error(501): Bad Latitude ... try again\n");
                exit(-1);
            }
            xlon = c2v(old_longitude2, 2);
            if (fabs(xlon + 999) < 0.001) {
                printf("Error(501): Bad Longitude ... try again\n");
                exit(-1);
            }



            char temp_latitude[50];
            char temp_longitude[50];
            char new_ellip_ht[15];

            strncpy(temp_latitude, "\0", sizeof (temp_latitude));
            strncpy(temp_longitude, "\0", sizeof (temp_longitude));
            strncpy(new_ellip_ht, "\0", sizeof (new_ellip_ht));
            char lat_dms[20];
            char lon_dms[20];
            char latdeg_c[12];
            char latmin_c[12];
            char latsec_c[12];
            char londeg_c[12];
            char lonmin_c[12];
            char lonsec_c[12];
            char cstddev[9];
            char cdist[9];

            int latdeg, latmin;
            int londeg, lonmin;
            double latsec, lonsec;

            // Initialize local variables
            strncpy(lat_dms, "\0", 20);
            strncpy(lon_dms, "\0", 20);
            strncpy(latdeg_c, "\0", 12);
            strncpy(latmin_c, "\0", 12);
            strncpy(latsec_c, "\0", 12);
            strncpy(londeg_c, "\0", 12);
            strncpy(lonmin_c, "\0", 12);
            strncpy(lonsec_c, "\0", 12);
            strncpy(cstddev, "\0", 9);
            strncpy(cdist, "\0", 9);
            dd_dms(xlat, lat_dms);

            //change to West longitude because HTDP only accepts West
            if (poseast == 1) { // 1 := true => is East lon)
                dd_dms((360.0 - xlon), lon_dms);
            } else {
                dd_dms(xlon, lon_dms);
            }

            strncpy(latdeg_c, &lat_dms[0], 3);
            strncpy(latmin_c, &lat_dms[3], 2);
            strncpy(latsec_c, &lat_dms[5], 8);

            strncpy(londeg_c, &lon_dms[0], 3);
            strncpy(lonmin_c, &lon_dms[3], 2);
            strncpy(lonsec_c, &lon_dms[5], 8);

            latdeg = atoi(latdeg_c);
            latmin = atoi(latmin_c);
            latsec = atof(latsec_c);

            londeg = atoi(londeg_c);
            lonmin = atoi(lonmin_c);
            lonsec = atof(lonsec_c);

            if ((latsec + 0.000001) >= 60.0) {
                latsec -= 60.0;
                ++latmin;
            }
            if (latmin >= 60) {
                latmin -= 60;
                ++latdeg;
            }

            if ((lonsec + 0.000001) >= 60.0) {
                lonsec -= 60.0;
                ++lonmin;
            }
            if (lonmin >= 60) {
                lonmin -= 60;
                ++londeg;
            }

            char dms_lat [16];
            char dms_lon [16];
            isHTDPUsed = 0;
            strncpy (dms_lat, "\0", sizeof(dms_lat));
            strncpy (dms_lon, "\0", sizeof(dms_lon));
            sprintf(dms_lat, "%3d %2d %8.5lf", latdeg, latmin, latsec);
            sprintf(dms_lon, "%3d %2d %8.5lf", londeg, lonmin, lonsec);
            //printf("isNAD83[%d] isPA11[%d] old_ellip_height[%s]\n",isNAD83,isPA11,old_ellip_height);
            if ((isNAD83 == 1 || isPA11 == 1 || isMA11 == 1) && (*old_ellip_height != '\0')) {
                //check if computing GEOID12B or others
                if (is_geoid12b){
                    //if input NAD83 compute same as geoid12b
                    if(strstr(inputType, "NAD83") != NULL) {
                        //printf ("input is in NAD83 do nothing\n");
                        strcpy (new_latitude[count], old_latitude);
                        strcpy (new_longitude[count], old_longitude);
                        new_ellip_height[count]  = atof(old_ellip_height);
                        //printf("NAD83 new_ellip_height[%lf]\n",new_ellip_height[count]);
                    }else{//else if input IGS08 convert IGS08 to 83 using HTDP
                            if (imodel >= 24 && imodel <= 25) {

                            }else{
                                HTDP_08_to_83(dms_lat, dms_lon, old_ellip_height, poseast, text, &new_latitude[count], &new_longitude[count], &new_ellip_ht);
                                isHTDPUsed = 1;
                                //printf("HTDP new_ellip_height[%lf]\n",new_ellip_height[count]);
                            }

                    }

                }else{
                    //if input NAD83 convert to NAD83 to IGS08 using HTDP
                    if (imodel >= 24 && imodel <= 25) {
                        strcpy (new_latitude[count], old_latitude);
                        strcpy (new_longitude[count], old_longitude);
                        new_ellip_height[count]  = atof(old_ellip_height);
                    }else{
                        if(strstr(inputType, "NAD83") != NULL) {
                            if (execHTDP == 0){
                                HTDP_83_to_08(dms_lat, dms_lon, old_ellip_height, poseast, text, &new_latitude[count], &new_longitude[count], &new_ellip_ht);
                                isHTDPUsed = 1;
                                execHTDP = 1;
                                //check if hemi='S' then add minus(-) sign
                                if (strstr(new_latitude[count], "S")){
                                    char dash[50];
                                    strncpy(dash,"\0", sizeof(dash));
                                    strcpy(dash,"-");
                                    strcat(dash,new_latitude[count]);
                                    strncpy(new_latitude[count],"\0", sizeof(new_latitude[count]));
                                    strncpy(new_latitude[count],dash,15);
                                }else{
                                    char dash[50];
                                    strncpy(dash,"\0", sizeof(dash));
                                    strcpy(dash,"");
                                    strcat(dash,new_latitude[count]);
                                    strncpy(new_latitude[count],"\0", sizeof(new_latitude[count]));
                                    strncpy(new_latitude[count],dash,14);
                                }

                                strcpy(saved_new_latitude, new_latitude[count]);
                                strcpy(saved_new_longitude, new_longitude[count]);
                                strcpy(saved_new_ellip_height, new_ellip_ht);
                            }else{
                                strcpy(new_latitude[count], saved_new_latitude);
                                strcpy(new_longitude[count], saved_new_longitude);
                                strcpy(new_ellip_ht, saved_new_ellip_height);
                            }

                        }else{//else if input IGS08 do nothing
                            //printf ("input is in IGS08 do nothing\n");
                            strcpy (new_latitude[count], old_latitude);
                            strcpy (new_longitude[count], old_longitude);
                            new_ellip_height[count]  = atof(old_ellip_height);
                        }
                    }

                }

                strcpy(temp_latitude, new_latitude[count]);
                strcpy (temp_longitude, new_longitude[count]);
                //printf("temp_latitude[%s]\n",temp_latitude);
                if (strcmp(temp_latitude, "") == 0){
                    geoid_height[count] = (double) -999.;
                }else{
                    xlat = c2v(temp_latitude, 1);
                    xlon = c2v(temp_longitude, 2);
                    // change back to East longitude because of HTDP
                    if (poseast  && isHTDPUsed) {
                        isEast = 1;
                        xlon = 360. - xlon;
                    }
                    new_ellip_height[count] = atof(new_ellip_ht);

                    //obtain ellipsoid height change by subtracting the original eht to the new eht
                    //dh = h_new - h_old ; h_new (igs08),  h_old (nad83)
                    dh[count] = atof(old_ellip_height) - new_ellip_height[count];  new_ellip_height[count] = atof(new_ellip_ht);

                    //obtain ellipsoid height change by subtracting the original eht to the new eht
                    //dh = h_new - h_old ; h_new (igs08),  h_old (nad83)
                    dh[count] = atof(old_ellip_height) - new_ellip_height[count];
                }
            }else{
                //check if computing GEOID12B or others
                if (is_geoid12b){
                    //if input IGS08 convert IGS08 to 83 using HTDP
                    if (imodel >= 24 && imodel <= 25) {

                    }else{
                        HTDP_08_to_83(dms_lat, dms_lon, old_ellip_height, poseast, text, &new_latitude[count], &new_longitude[count], &new_ellip_ht);
                        isHTDPUsed = 1;
                        //printf("HTDP new_ellip_ht[%s]\n",new_ellip_ht);
                    }
                }else {
                    if (imodel >= 24 && imodel <= 25) {

                    }else{
                        if (execHTDP == 0) {
                           HTDP_08_to_08(dms_lat, dms_lon, old_ellip_height, poseast, text, &new_latitude[count], &new_longitude[count], &new_ellip_ht);
                           isHTDPUsed = 1;
                           execHTDP = 1;

                           //remove hemisphere 'N' from HTDP output
                           char dash[50];
                           strncpy(dash,"\0", sizeof(dash));
                           strcpy(dash,"");
                           strcat(dash,new_latitude[count]);
                           strncpy(new_latitude[count],"\0", sizeof(new_latitude[count]));
                           strncpy(new_latitude[count],dash,14);
                           strcpy(saved_new_latitude, new_latitude[count]);
                           strcpy(saved_new_longitude, new_longitude[count]);
                           strcpy(saved_new_ellip_height, new_ellip_ht);
                       }else{
                               strcpy(new_latitude[count], saved_new_latitude);
                               strcpy(new_longitude[count], saved_new_longitude);
                               strcpy(new_ellip_ht, saved_new_ellip_height);
                       }
                    }
                }
                strcpy(temp_latitude, new_latitude[count]);
                strcpy(temp_longitude, new_longitude[count]);
                //printf("temp_latitude[%s]\n",temp_latitude);
                if (strcmp(temp_latitude, "") == 0) {
                    geoid_height[count] = (double) - 999.;
                } else {
                    xlat = c2v(temp_latitude, 1);
                    xlon = c2v(temp_longitude, 2);
                    // change back to East longitude because of HTDP
                    if (poseast && isHTDPUsed) {
                        isEast = 1;
                        xlon = 360. - xlon;
                    }
                    new_ellip_height[count] = atof(new_ellip_ht);

                    //obtain ellipsoid height change by subtracting the original eht to the new eht
                    //dh = h_new - h_old ; h_new (igs08),  h_old (nad83)
                    dh[count] = atof(old_ellip_height) - new_ellip_height[count];
                    new_ellip_height[count] = atof(new_ellip_ht);

                    //obtain ellipsoid height change by subtracting the original eht to the new eht
                    //dh = h_new - h_old ; h_new (igs08),  h_old (nad83)
                    dh[count] = atof(old_ellip_height) - new_ellip_height[count];
                }
            }

            // Force input West longitude to be positive East
            if (poseast == 0) {       // C std: 0 = false; is NOT east lon
                xlon = 360. - xlon;
            }

            // load lat/lon struct to 1-D array
            DATASET1 thisSet;
            strncpy(thisSet.text, "\0", 50);
            thisSet.lat     = xlat;
            thisSet.lon     = xlon;
            thisSet.poseast = poseast;
            strcpy(thisSet.text, text);
            vec_data[kount] = thisSet;
            ++kount;
            if (kount >= mem_limit) {
                mem_limit += MEM_STEP;
                vec_data =
                    (DATASET1*)realloc(vec_data, mem_limit * sizeof(DATASET1));
                if (vec_data == NULL) {
                    printf("Out of system memory - allocation fails\n");
                    exit(-1);
                }
            }

        //} while (jj == 0);

    }//~(iinput == 1) keyboard prompt

    // ---------------------------------------------
    // Now have the lat/lon pair(s) in a data vector
    // ---------------------------------------------


    // ===== COMMON TO NON-BLUEBOOK INPUT TYPES ============

    // ===== INTERPOLATE GEOID VALUE =======================

    // iterate thru the array -----
    for (ii = 0; ii < kount; ++ii) {
        xlat = vec_data[ii].lat;
        xlon = vec_data[ii].lon;

        if (DEBUG != 0) {
            printf("kount = %d  ii = %d    xlat  = %lf   xlon = %lf\n",
                    kount, ii, xlat, xlon);
        }

        // If the lat/lon values came back as -999, set the
        // geoid value to -999 and skip the interpolation
        if (xlat == -999. || xlon == -999.) {
            geoidHt = (double) -999.;
            stddev = (double) -999.;
            distance = (double) -999.;
            continue;
        } else {

            // Find which geoid file to use, based on the lat/lon input
            kk = which1( xlat, xlon, nfiles, kk, imodel, vec_fnames, vec_hdr, vec_ifp, count );
/*
            if (kk == 3){
                isPR = 1;
            }else{
                isPR = 0;
            }
*/

            if (DEBUG != 0) { printf("kk = %d \n", kk); }

            // If the point isn't in any of our grid areas, set to -999
            if (kk == -1) {
                geoidHt = (double) -999.;
                stddev = (double) -999.;
                distance = (double) -999.;
                geoid_height[count] = geoidHt;
                ortho_height[count] = 0.0;

            // Otherwise, do the interpolation
            } else {
                if (DEBUG != 0) {
                    printf("xlat = %lf  xlon = %lf  kk = %d\n ", xlat,xlon,kk);
                }

                geoidHt = interg(xlat, xlon, vec_hdr, vec_ifp, kk, count );

                //scale the geoid height into the new ellipsoid datum from NAD83 to IGS08, or vice versa
                if ((isNAD83 || isPA11 || isMA11) && is08_to_83){
                    //N_new = N_old - dh
                    //printf("geoid[%lf] - dh[%lf]\n",geoidHt,dh[count]);
                    geoidHt = geoidHt - dh[count];
                }
                geoid_height[count] = geoidHt;


                //compute orthometric height if ellipsoid height exist
                //H = h - N ; ortho_ht = ellip_ht - geoid_ht
                //check if computing GEOID12B
                if (is_geoid12b){
                    if (*old_ellip_height != '\0') {
                        if(strstr(inputType, "NAD83") != NULL) {
                            //if NAD83 then just use input ellip_height
                            ortho_height[count] = atof(old_ellip_height) - geoidHt;
                        }else{
                            //use converted ellip_height
                            //printf("eh[%lf] - geoid[%lf]\n",new_ellip_height[0],geoidHt);
                            ortho_height[count] = new_ellip_height[count] - geoidHt;
                        }
                    }
                }else{
                    if (*old_ellip_height != '\0') {
                        if(strstr(inputType, "NAD83") != NULL) {
                            //use converted ellip_height
                            ortho_height[count] = new_ellip_height[count] - geoidHt;
                        }else{
                            //if IGS08 then just use input ellip_height
                            ortho_height[count] = atof(old_ellip_height) - geoidHt;
                        }
                    }
                }

                if (DEBUG != 0) { printf("geoidHt = %lf \n", geoidHt ); }

                //  If we have variance and distance grids available,
                //  interpolate from those as well
                stddev = (double) -999.;
                distance = (double) -999.;


                if ( nvdfiles > 0 ) {
                   //printf("xlat[%lf] xlon[%lf]\n",xlat,xlon);
                   mm = which1( xlat, xlon, nvdfiles, mm, imodel, vdar_fnames,
                                vdar_hdr, vdar_ifp, dummyCount );
                   if ( mm != -1 ) {
                      //printf("vdariance = interg_idw\n");
                      vdariance = interg_idw(xlat, xlon, vdar_hdr, vdar_ifp, mm );
                      DN = vdariance;
                      if (DEBUG != 0) { printf("DN = %lf \n", DN ); }
                   }
                }

                if ( nvfiles > 0 ) {
                   //printf("xlat[%lf] xlon[%lf]\n",xlat,xlon);
                   mm = which1( xlat, xlon, nvfiles, mm, imodel, var_fnames,
                                var_hdr, var_ifp, dummyCount );
                   if ( mm != -1 ) {
                      //printf("variance = interg_idw\n");
                      variance = interg_idw(xlat, xlon, var_hdr, var_ifp, mm );
                      STD = 1.96 * variance;
                      if (DEBUG != 0) { printf("STD = %lf \n", STD ); }
                   }
                }

                //printf("SAGout: %lf +/-%lf m %lf \n",geoidHt,stddev,distance);

                ++keep;
                ave += geoidHt;
                rms += geoidHt * geoidHt;
                if (xlat > xn)  xn = xlat;
                if (xlat < xs)  xs = xlat;
                if (xlon > xe)  xe = xlon;
                if (xlon < xw)  xw = xlon;

                if (geoidHt < minght) {
                    minght = geoidHt;
                    minlat = xlat;
                    minlon = xlon;
                }
                if (geoidHt > maxght) {
                    maxght = geoidHt;
                    maxlat = xlat;
                    maxlon = xlon;
                }
            }

            if (strcmp(ofyn, "Y") == 0 || strcmp(ofyn, "y") == 0) {
                if      (iform == 1)  ff1out(ofp, vec_data[ii], geoidHt, imodel, stddev, distance);
                else if (iform == 2)  ff2out(ofp, vec_data[ii], geoidHt, imodel, stddev, distance);
                // else if (iform == 3)  // bluebook - do nothing -----

                else if (iinput == 1) ff1out(ofp, vec_data[ii], geoidHt, imodel, stddev, distance);

            }//~if(ofyn)

            if( (iinput == 1) && (strcmp(ofyn, "Y") != 0) ){
                ff4out(ofp, vec_data[ii], geoidHt, imodel, stddev, distance);
            }


            // if( (iinput == 1) && (outfil == 1) )
            //     ff4out(ofp, vec_data[ii], geoidHt, imodel, stddev, distance);

        }//~if(xlat == -999. || xlon == -999.)

    }//~for(ii)

    free( (void*)vec_data );

    //printf("DONE\n");
    // Finally, write out the record to screen and possibly to an output file

    // ===== OUTPUT =======================
    // write output to screen and possibly to an output file
    // go get another input dataset


    // =============================================================
    // This is for input file format = 1
    // When finished, give a little report to the screen and end program.
    //
    ave = ave/keep;
    rms = sqrt(rms/keep);
    if (keep > 1) {
        fact = ((double)(keep)/(double)(keep-1));
        std  = sqrt( fact*(rms*rms - ave*ave) );
    } else {
        std = 0;
    }

/*
    printf("\n%s\n", dash70);
    printf("\
FINAL REPORT: \n\
  Number of points input: %8d \n\
  Number of good points : %8d \n\
  Northernmost Latitude : %10.6lf  \n\
  Southernmost Latitude : %10.6lf \n", kount, keep, xn, xs);

    if (poseast == 0) { // poseast = 0  := Westlon
        printf("\
  Westernmost Longitude : %10.6lf \n\
  Easternmost Longitude : %10.6lf \n\
  Minimum Geoid Height  :  %8.3lf \n\
    Lat/Lon of Minimum  : %10.6lf  %10.6lf \n\
  Maximum Geoid Height  :  %8.3lf \n\
    Lat/Lon of Maximum  : %10.6lf  %10.6lf \n\
  Average Geoid Height  :  %8.3lf \n\
  Standard Deviation    :  %8.3lf \n\
  Root Mean Square      :  %8.3lf \n",
            (360. - xw), (360. - xe), minght, minlat, (360. - minlon),
            maxght, maxlat, (360. - maxlon), ave, std, rms);

    } else {       // using East longitude

        printf("\
  Westernmost Longitude : %10.6lf \n\
  Easternmost Longitude : %10.6lf \n\
  Minimum Geoid Height  :  %8.3lf \n\
    Lat/Lon of Minimum  : %10.6lf  %10.6lf \n\
  Maximum Geoid Height  :  %8.3lf \n\
    Lat/Lon of Maximum  : %10.6lf  %10.6lf \n\
  Average Geoid Height  :  %8.3lf \n\
  Standard Deviation    :  %8.3lf \n\
  Root Mean Square      :  %8.3lf \n",
            xw, xe, minght, minlat, minlon,
            maxght, maxlat, maxlon, ave, std, rms);

    }
*/

        //printf("CLOSE vec_ifp\n");
        int i = 0;

        for (i = 0; i < 25; ++i) {
			//printf("i[%d]\n",i);
            if (vec_ifp[count][i]!= NULL){
               //printf("vec_ifp[%s] count[%d] i[%d]\n",vec_ifp[count][i], count, i);
               fclose(vec_ifp[count][i]);
               vec_ifp[count][i] = NULL;
            }
        }


        //printf("CLOSE vdar_ifp\n");
        //printf("imodel[%d]\n",imodel);
        if (imodel >= 24 && imodel <= 27) {
            for (i = 0; i < 3; ++i){
               if (vdar_ifp[i] != NULL){
                  //printf("vdar_ifp[i][%s]\n",vdar_ifp[i]);
                  fclose(vdar_ifp[i]);
			  }
            }
        }

        //printf("CLOSE var_ifp\n");
        if (imodel == 27) {
            for (i = 0; i < 1; ++i){
               if (var_ifp[i] != NULL)
                  fclose(var_ifp[i]);
            }
        }

        //printf("CLOSE ifp\n");
    if (iinput == 2) {  // input by file
        fclose(ifp);
    }
        //printf("CLOSE ofp\n");
    if (strcmp(ofyn, "Y") == 0 || strcmp(ofyn, "y") == 0) {
        fclose(ofp);
    }

    return(0);
}

/**
 * call HTDP to transform NAD83 to IGS08
 *
 * @param latitude lat from user input
 * @param longitude lon from user input
 * @param ellip_ht ellipsoid height from user input
 * @param poseast is true if position is east longitude
 * @param text station name
 * @param new_latitude converted lat being returned
 * @param new_longitude converted lon being returned
 * @param new_ellip_ht convert ellipsoid height being returned
 * @return
 */
int HTDP_83_to_08 (char *latitude, char *longitude, char *ellip_ht, int poseast, char *text, char *new_latitude, char *new_longitude, char *new_ellip_ht){
    FILE *ifp;
    FILE *ofp;
    char ifile[256];
    char ofile[256];
    char datetime[256];
    char station_name[256];
    char htdp_in[256];
    char htdp_out[256];
    char line[128];
    char dummy[2];
    int datum=1; //1 - NAD83(2011), 2 - NAD83(PA11), 3 - NAD83 (MA11)

    if (isNAD83){
        datum = 1;
    }if (isPA11){
        datum = 2;
    }else {
        datum = 3;
    }

    //create input file for HTDP using process id
    strncpy(datetime, "\0", sizeof (datetime));
    strncpy(ifile, "\0", sizeof (ifile));
    strncpy(ofile, "\0", sizeof (ofile));
    strncpy(dummy, "\0", sizeof (dummy));
    strncpy(htdp_in, "\0", sizeof (htdp_in));
    strncpy(htdp_out, "\0", sizeof (htdp_out));
    strcpy (dummy, ".");
    strcpy (station_name, text);

#ifdef NGS_PC_ENV
	strcpy(htdp_in, "htdp_in");
	strcpy(htdp_out, "htdp_out");
	strcpy(ifile, htdp_in);
	strcpy(ofile, htdp_out);
#else
	strcpy(htdp_in, "/tmp/htdp_in");
	strcpy(htdp_out, "/tmp/htdp_out");
	strcpy(ifile, htdp_in);
	strcpy(ofile, htdp_out);
	cat_date(dummy, datetime);
	strcat (ifile, datetime);
        strcat (ofile, datetime);
#endif

    if ((ifp = fopen(ifile, "w")) == NULL) {
        printf("\n ABORT: Can not open %s.\n", ifile);
#ifndef NGS_PC_ENV
        exit(-1);
#else
        exit(1);
#endif
    }
    fprintf(ifp,"\n");
    fprintf(ifp,"4\n");
    fprintf(ifp,"%s\n", ofile);
    fprintf(ifp,"%d\n",datum);
    fprintf(ifp,"24\n");
    fprintf(ifp,"2\n");
    fprintf(ifp,"2010.0\n");
    fprintf(ifp,"2\n");
    fprintf(ifp,"2020.0\n");
    fprintf(ifp,"1\n");
    fprintf(ifp,"%s\n", station_name);
    fprintf(ifp,"1\n");
    fprintf(ifp,"%s\n", latitude);
    fprintf(ifp,"%s\n", longitude);
    fprintf(ifp,"%s\n", ellip_ht);
    fprintf(ifp,"0\n");
    fprintf(ifp,"n\n");
    fprintf(ifp,"0\n");

    fflush(ifp);
    fclose(ifp);

    char htdp_cmd[75];
    strncpy(htdp_cmd, "\0", 75);

#ifdef NGS_PC_ENV
	strcpy(htdp_cmd, pcDir);
	strcpy(htdp_cmd, "htdp.exe < ");
	strcat(htdp_cmd, ifile);
    strcat(htdp_cmd, " > NUL");
    system(htdp_cmd)
#else
	strcpy(htdp_cmd, "htdp < ");
	strcat(htdp_cmd, ifile);
    strcat(htdp_cmd, " > /dev/null");
    system(htdp_cmd)
#endif

    ;

    if ((ofp = fopen(ofile, "r")) == NULL) {
        printf("\n ABORT: Can not open %s\n", ofile);
#ifndef NGS_PC_ENV
        exit(-1);
#else
        exit(1);
#endif
    }

    strncpy(line, "\0", sizeof(line));
    strncpy (new_ellip_ht, "\0", sizeof(new_ellip_ht));
    while ((int) fgets(line, 127, ofp) != 0) {
        //printf ("out - %s", line);
        if(strstr(line, "LATITUDE") != NULL) {
            strncpy(new_latitude, line + 34, 18);
            trim_c( new_latitude, 'b');
        }
        if(strstr(line, "LONGITUDE") != NULL) {
            strncpy(new_longitude, line + 32, 18);
            trim_c( new_longitude, 'b');
        }
        if(strstr(line, "ELLIP") != NULL) {
            strncpy(new_ellip_ht, line + 38, 14);
            trim_c( new_ellip_ht, 'b');
        }
    }
    if(new_ellip_ht[0] == '\0') {
        printf("HTDP conversion error - check your inputs\n");
        exit(1);
    }
    fclose(ofp);
    remove(ifile);
    remove(ofile);
    return (1);
}

/**
 * call HTDP to transform IGS08 epoch 2005.00 to IGS08 epoch 2022.00
 *
 * @param latitude lat from user input
 * @param longitude lon from user input
 * @param ellip_ht ellipsoid height from user input
 * @param poseast is true if position is east longitude
 * @param text station name
 * @param new_latitude converted lat being returned
 * @param new_longitude converted lon being returned
 * @param new_ellip_ht convert ellipsoid height being returned
 * @return
 */
int HTDP_08_to_08 (char *latitude, char *longitude, char *ellip_ht, int poseast, char *text, char *new_latitude, char *new_longitude, char *new_ellip_ht){
    FILE *ifp;
    FILE *ofp;
    char ifile[256];
    char ofile[256];
    char datetime[256];
    char station_name[256];
    char htdp_in[256];
    char htdp_out[256];
    char line[128];
    char dummy[2];

    //create input file for HTDP using process id
    strncpy(datetime, "\0", sizeof (datetime));
    strncpy(ifile, "\0", sizeof (ifile));
    strncpy(ofile, "\0", sizeof (ofile));
    strncpy(dummy, "\0", sizeof (dummy));
    strncpy(htdp_in, "\0", sizeof (htdp_in));
    strncpy(htdp_out, "\0", sizeof (htdp_out));
    strcpy (dummy, ".");
    strcpy (station_name, text);

#ifdef NGS_PC_ENV
	strcpy(htdp_in, "htdp_in");
	strcpy(htdp_out, "htdp_out");
	strcpy(ifile, htdp_in);
	strcpy(ofile, htdp_out);
#else
	strcpy(htdp_in, "/tmp/htdp_in");
	strcpy(htdp_out, "/tmp/htdp_out");
	strcpy(ifile, htdp_in);
	strcpy(ofile, htdp_out);
	cat_date(dummy, datetime);
	strcat (ifile, datetime);
        strcat (ofile, datetime);
#endif

    if ((ifp = fopen(ifile, "w")) == NULL) {
        printf("\n ABORT: Can not open %s.\n", ifile);
#ifndef NGS_PC_ENV
        exit(-1);
#else
        exit(1);
#endif
    }
    fprintf(ifp,"\n");
    fprintf(ifp,"4\n");
    fprintf(ifp,"%s\n", ofile);
    fprintf(ifp,"24\n");
    fprintf(ifp,"24\n");
    fprintf(ifp,"2\n");
    fprintf(ifp,"2010.0\n");
    fprintf(ifp,"2\n");
    fprintf(ifp,"2020.0\n");
    fprintf(ifp,"1\n");
    fprintf(ifp,"%s\n", station_name);
    fprintf(ifp,"1\n");
    fprintf(ifp,"%s\n", latitude);
    fprintf(ifp,"%s\n", longitude);
    fprintf(ifp,"%s\n", ellip_ht);
    fprintf(ifp,"0\n");
    fprintf(ifp,"n\n");
    fprintf(ifp,"0\n");

    fflush(ifp);
    fclose(ifp);

    char htdp_cmd[75];
    strncpy(htdp_cmd, "\0", 75);

#ifdef NGS_PC_ENV
	strcpy(htdp_cmd, pcDir);
	strcpy(htdp_cmd, "htdp.exe < ");
	strcat(htdp_cmd, ifile);
    strcat(htdp_cmd, " > NUL");
    system(htdp_cmd)
#else
	strcpy(htdp_cmd, "htdp < ");
	strcat(htdp_cmd, ifile);
    strcat(htdp_cmd, " > /dev/null");
    system(htdp_cmd)
#endif

    ;

    if ((ofp = fopen(ofile, "r")) == NULL) {
        printf("\n ABORT: Can not open %s\n", ofile);
#ifndef NGS_PC_ENV
        exit(-1);
#else
        exit(1);
#endif
    }

    strncpy(line, "\0", sizeof(line));
    strncpy (new_ellip_ht, "\0", sizeof(new_ellip_ht));
    while ((int) fgets(line, 127, ofp) != 0) {
        //printf ("out - %s", line);
        if(strstr(line, "LATITUDE") != NULL) {
            strncpy(new_latitude, line + 34, 18);
            trim_c( new_latitude, 'b');
        }
        if(strstr(line, "LONGITUDE") != NULL) {
            strncpy(new_longitude, line + 32, 18);
            trim_c( new_longitude, 'b');
        }
        if(strstr(line, "ELLIP") != NULL) {
            strncpy(new_ellip_ht, line + 38, 14);
            trim_c( new_ellip_ht, 'b');
        }
    }
    if(new_ellip_ht[0] == '\0') {
        printf("HTDP conversion error - check your inputs\n");
        exit(1);
    }
    fclose(ofp);
    remove(ifile);
    remove(ofile);
    return (1);
}

/**
 * call HTDP to transform IGS08 to NAD83
 *
 * @param latitude lat from user input
 * @param longitude lon from user input
 * @param ellip_ht ellipsoid height from user input
 * @param poseast is true if position is east longitude
 * @param text station name
 * @param new_latitude converted lat being returned
 * @param new_longitude converted lon being returned
 * @param new_ellip_ht convert ellipsoid height being returned
 * @return
 */
int HTDP_08_to_83 (char *latitude, char *longitude, char *ellip_ht, int poseast, char *text, char *new_latitude, char *new_longitude, char *new_ellip_ht){
    FILE *ifp;
    FILE *ofp;
    char ifile[256];
    char ofile[256];
    char datetime[256];
    char station_name[256];
    char htdp_in[256];
    char htdp_out[256];
    char line[128];
    char dummy[2];
    int datum=1; //1 - NAD83(2011), 2 - NAD83(PA11), 3 - NAD83(MA11)

    //check if coordinate is on pacific plate - NAD83 (PA11)
    //else it's NAD83 (2011)
    char dms_lat [16];
    char dms_lon [16];
    strcpy(dms_lat, latitude);
    strcpy(dms_lon, longitude);
    if (am_I_on_the_Pacific_plate(dms_dd (dms_lat), dms_dd (dms_lon))){
        //printf ("YESSSSSS\n");
        datum = 2;
    }else{
        //printf ("NOOOOOOO\n");
        datum = 1;
    }

    //create input file for HTDP using process id
    strncpy(datetime, "\0", sizeof (datetime));
    strncpy(ifile, "\0", sizeof (ifile));
    strncpy(ofile, "\0", sizeof (ofile));
    strncpy(dummy, "\0", sizeof (dummy));
    strncpy(htdp_in, "\0", sizeof (htdp_in));
    strncpy(htdp_out, "\0", sizeof (htdp_out));
    strcpy (dummy, ".");
    strcpy (station_name, text);

#ifdef NGS_PC_ENV
	strcpy(htdp_in, "htdp_in");
	strcpy(htdp_out, "htdp_out");
	strcpy(ifile, htdp_in);
	strcpy(ofile, htdp_out);
#else
	strcpy(htdp_in, "/tmp/htdp_in");
	strcpy(htdp_out, "/tmp/htdp_out");
	strcpy(ifile, htdp_in);
	strcpy(ofile, htdp_out);
	cat_date(dummy, datetime);
	strcat (ifile, datetime);
    strcat (ofile, datetime);
#endif


    if ((ifp = fopen(ifile, "w")) == NULL) {
        printf("\n ABORT: Can not open %s\n", ifile);
#ifndef NGS_PC_ENV
        exit(-1);
#else
        exit(1);
#endif
    }
    fprintf(ifp,"\n");
    fprintf(ifp,"4\n");
    fprintf(ifp,"%s\n", ofile);
    fprintf(ifp,"24\n");
    fprintf(ifp,"%d\n", datum);
    fprintf(ifp,"2\n");
    fprintf(ifp,"2020.0\n");
    fprintf(ifp,"2\n");
    fprintf(ifp,"2010.0\n");
    fprintf(ifp,"1\n");
    fprintf(ifp,"%s\n", station_name);
    fprintf(ifp,"1\n");
    fprintf(ifp,"%s\n", latitude);
    fprintf(ifp,"%s\n", longitude);
    fprintf(ifp,"%s\n", ellip_ht);
    fprintf(ifp,"0\n");
    fprintf(ifp,"n\n");
    fprintf(ifp,"0\n");

    fflush(ifp);
    fclose(ifp);

    char htdp_cmd[75];
    strncpy(htdp_cmd, "\0", 75);

#ifdef NGS_PC_ENV
	strcpy(htdp_cmd, pcDir);
	strcat(htdp_cmd, "htdp.exe < ");
	strcat(htdp_cmd, ifile);
	strcat(htdp_cmd, " > NUL");
    system(htdp_cmd);
#else
	strcpy(htdp_cmd, "htdp < ");
	strcat(htdp_cmd, ifile);
	strcat(htdp_cmd, " > /dev/null");
    system(htdp_cmd);
#endif


    if ((ofp = fopen(ofile, "r")) == NULL) {
        printf("\n ABORT: Can not open %s\n", ofile);
#ifndef NGS_PC_ENV
        exit(-1);
#else
        exit(1);
#endif
    }

    strncpy(line, "\0", sizeof(line));
    strncpy (new_ellip_ht, "\0", sizeof(new_ellip_ht));
    while ((int) fgets(line, 127, ofp) != 0) {
        //printf ("out - %s", line);
        if(strstr(line, "LATITUDE") != NULL) {
            strncpy(new_latitude, line + 32, 18);
            trim_c( new_latitude, 'b');
        }
        if(strstr(line, "LONGITUDE") != NULL) {
            strncpy(new_longitude, line + 32, 18);
            trim_c( new_longitude, 'b');
        }
        if(strstr(line, "ELLIP") != NULL) {
            strncpy(new_ellip_ht, line + 38, 14);
            trim_c( new_ellip_ht, 'b');
        }
    }
    if(new_ellip_ht[0] == '\0') {
        printf("HTDP conversion error - check your inputs\n");
        exit(1);
    }
    fclose(ofp);
    remove(ifile);
    remove(ofile);
    return (1);
}

char * strtok_single (char * str, char const * delims)
{
  static char  * src = NULL;
  char  *  p,  * ret = 0;

  if (str != NULL)
    src = str;

  if (src == NULL)
    return NULL;

  if ((p = strpbrk (src, delims)) != NULL) {
    *p  = 0;
    ret = src;
    src = ++p;

  } else if (*src) {
    ret = src;
    src = NULL;
  }

  return ret;
}
