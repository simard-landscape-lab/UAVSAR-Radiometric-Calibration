#include <iostream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <fstream>
#include <complex>
#include <time.h>
#include <float.h>
#include <cstring>
#include "optionparser.h"
#include "load_ann.h"

using namespace std;

extern char *optarg;
extern int optopt;


string SplitFilename (const string& str)
{
    size_t found;
    string path;

    found=str.find_last_of("/\\");
    
    if (found < str.size())
        path = str.substr(0,found+1);
    else
        path = "";

    return path;
}


struct Arg: public option::Arg
{
    static void printError(const char* msg1, const option::Option& opt, const char* msg2)
    {
        fprintf(stderr, "%s", msg1);
        fwrite(opt.name, opt.namelen, 1, stderr);
        fprintf(stderr, "%s", msg2);
    }
    
    static option::ArgStatus Unknown(const option::Option& option, bool msg)
    {
        if (msg) printError("Unknown option '", option, "'\n");
        return option::ARG_ILLEGAL;
    }
    
    static option::ArgStatus Required(const option::Option& option, bool msg)
    {
        if (option.arg != 0)
            return option::ARG_OK;
        
        if (msg) printError("Option '", option, "' requires an argument\n");
        return option::ARG_ILLEGAL;
    }
    
    static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
    {
        if (option.arg != 0 && option.arg[0] != 0)
            return option::ARG_OK;
        
        if (msg) printError("Option '", option, "' requires a non-empty argument\n");
        return option::ARG_ILLEGAL;
    }
    
    static option::ArgStatus Numeric(const option::Option& option, bool msg)
    {
        char* endptr = 0;
        if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
        if (endptr != option.arg && *endptr == 0)
            return option::ARG_OK;
        
        if (msg) printError("Option '", option, "' requires a numeric argument\n");
        return option::ARG_ILLEGAL;
    }
};


enum  optionIndex { UNKNOWN, HELP, OUT, CORR, AREA, TRANSIN, TRANSOUT, SIM, LOOK, SLOPE, MASK, RATIO };
const option::Descriptor usage[] =
{
    {UNKNOWN, 0, "", "",Arg::None, "Usage: uavsar_calib [-c vegetation_lut] [-a output_area] [-t trans_in] [-u trans_out] [-i local_incidence_out] [-l look_angle_out] [-s slope_angle_out] [-m mask_out] <ann file> <pol> <output intensity image>\n\n"
        "Required Arguments:" },
    {UNKNOWN, 0, "", "",Arg::None, "  <ann file>\tAnnotation file.\n  <pol>\t4-letter polarization string (HHHH, HVHV, or VVVV).\n  <output file>\tDestination filename for radiometrically calibrated intensity image.\n\n"
        "Optional Arguments:" },
    {HELP, 0,"h", "help",Arg::None,"  -h  \tPrint usage and exit." },
    {OUT, 0,"o", "out",Arg::Required, "  -o <output file>  \tOptional flag to save corrected intensity image." },
    {CORR, 0,"c", "corr",Arg::Required, "  -c <lut file>  \tOptional flag to perform LUT vegetation correction using provided LUT file." },
    {AREA, 0,"a", "area",Arg::Required, "  -a <area file>  \tOptional flag to save illuminated area image in RDC coordinates." },
    {TRANSIN, 0,"t", "tin",Arg::Required, "  -t <input transformation lut>  \tOptional flag to specify input transformation look up table." },
    {TRANSOUT, 0,"u", "tout",Arg::Required, "  -u <output transformation lut>  \tOptional flag to save output transformation look up table (for geocoding)." },
    {SIM, 0,"i","sim",Arg::Required, "  -i <local incidence file>  \tOptional flag to save local incidence angle map." },
    {LOOK, 0,"l","look",Arg::Required, "  -l <look file>  \tOptional flag to save look angle map." },
    {SLOPE, 0,"s","slope",Arg::Required, "  -s <slope file>  \tOptional flag to save output range-facing terrain slope angle map." },
    {MASK, 0, "m", "mask",Arg::Required, "  -m <mask file>  \tOptional flag to create a validity mask file which shows pixels where the radiometric calibration could not be performed.  Only valid for vegetation LUT correction using the -c option." },
    {RATIO, 0, "r", "ratio",Arg::Required, "  -r <ratio file>  \tOptional flag to create a ratio file which contains the ratio between the calibrated and uncalibrated images." },
    {UNKNOWN, 0, "", "",Arg::None, "\nExample Usage:\n"
        "  uavsar_calib -c caltbl_NewHampshire_WhiteMountain_HH.flt -u geomap.trans Brtlet_07101_09061_001_090814_L090_CX_01.ann HHHH Brtlet_HHHH_Cal.mlc "},
    {0,0,0,0,0,0}
};



int main(int argc, char* argv[]){
    argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
    option::Stats  stats(usage, argc, argv);
    std::vector<option::Option> options(stats.options_max);
    std::vector<option::Option> buffer(stats.buffer_max);
    option::Parser parse(usage, argc, argv, &options[0], &buffer[0]);
    
    if (parse.error())
        return 1;
    
    if (options[HELP] || argc == 0) {
        option::printUsage(std::cout, usage);
        return 0;
    }
    

    int iter, max_iter = 30, ix1, ix2, iy1, iy2, area_flag = 0, LUTin_flag = 0, LUTout_flag = 0, sim_flag = 0, correct_flag = 0, look_flag = 0, slope_flag = 0, mask_flag = 0, ratio_flag = 0,
      error_flag = 0, poly_method, cos_flag = 0,pol=5,e_look,e_slope,size;

    float Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Zavg, azpix, ranpix, p, q, deltaDEM_lat, deltaDEM_lon, xbound, ybound, x1, x2, y1, y2, cs, ss, tempout, h, r_area_fe, dist, fx1, fy1;
    
    float min_correction_ratio = 0.001, max_correction_ratio = 1000.0, void_correction_val = -1.0;
        
    double ta, tcen, satdist, earth_radius, sat_alt, vsat, velx, alpha, i_cur, j_cur, lat, lon, lvm,  
      slope, aspect, area, area_ref, temp, inc_cor, inc_tol, inc_ltol,slt_range, r_x1, r_x2, slope_r, slope_a, r_l3, slope_actual_r, slope_actual_a, r_look, theta_c, antcor;

    XYZ satxyz, satuvw, Xpix, lkv, SCH, SCH2, temp_raU, nI, nE, look_sch, nL;

    par_struct par;
    peg_struct peg;

    vector<float> LUTcpx(2,0);

    string area_out, LUT_flout, LUT_flin, sim_name, name_orbit, veg_in, look_name, slope_name, mask_name, diff_name;

    ifstream LUTin, VegTablefile;
    ofstream LUTout, areaRDCout, ampout, sim_flout, look_out, slope_out, mask_out, ratio_out;

    char *temp_char = NULL;
    
    time_t timerstart, timerend;    

    //Start timer
    timerstart = time (NULL);

    //-----------------------------------   Get command line arguments   ---------------------------------------
    if (argc < 3)
      error_flag = 1;
    else {
      par.ann = argv[argc-3];
      par.pol = argv[argc-2];

      if (!strcmp(argv[argc-2],"HHHH")) {
        cout <<"Calibrating for polarization HHHH"<<endl;
        pol = 1;
      }

      if (!strcmp(argv[argc-2],"HVHV")) {
        cout <<"Calibrating for polarization HVHV"<<endl;
        pol = 2;
      }

      if (!strcmp(argv[argc-2],"VVVV")) {
        cout <<"Calibrating for polarization VVVV"<<endl;
        pol = 3;
      }
        
      //Filename for correct intensity image.
      correct_flag = 1;
      par.cor_out = argv[argc-1];


      //Get optional command line arguments
      for (int i = 0; i < parse.optionsCount(); ++i) {
          option::Option& opt = buffer[i];
          switch (opt.index()) {
              case HELP:
                  // not possible, because handled further above and exits the program
                  break;
              case CORR:
                  cos_flag = 1;
                  veg_in = opt.arg;
                  break;
              case AREA:
                  area_flag = 1;
                  area_out = opt.arg;
                  break;
              case TRANSIN:
                  LUTin_flag = 1;
                  LUT_flin = opt.arg;
                  break;
              case TRANSOUT:
                  LUTout_flag = 1;
                  LUT_flout = opt.arg;
                  break;
              case SIM:
                  sim_flag = 1;
                  sim_name = opt.arg;
                  break;
              case LOOK:
                  look_flag = 1;
                  look_name = opt.arg;
                  break;
              case SLOPE:
                  slope_flag = 1;
                  slope_name = opt.arg;
                  break;
              case MASK:
                  mask_flag = 1;
                  mask_name = opt.arg;
                  break;
              case RATIO:
                  ratio_flag = 1;
                  diff_name = opt.arg;
                  break;            
          }
      }
        
    }
    
    if (error_flag){
      option::printUsage(std::cout, usage);
      return 0;
    }

    cout << "\nUAVSAR radiometric calibration software designed and written by Marc Simard and Bryan V. Riel.\n\n";
    cout << "\nCopyright 2010, by the California Institute of Technology. ALL RIGHTS RESERVED. \n";
    cout << "\nUnited States Government Sponsorship acknowledged. \n";
    cout << "\nAny commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.  This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons.\n\n";

    cout <<"\n-----------------------------------------------------------------------\n";


    cout << "\nPerforming facet model correction......\n\n";

    //----------------------------------------   Pre-processing steps  -------------------------------------------
    
    //Load information from annotation file
    load_ann(par, peg);
    xbound = (float)par.width-1; //Bounds for valid RDC coordinates
    ybound = (float)par.height-1;
    
    // Get path to annotation file.  Assume MLC, HGT, etc. files are in same folder as ANN.
    string path = SplitFilename(par.ann.c_str());

    //Create input and output file identifiers
    string mlcfile = path + par.mlc;
    string hgtfile = path + par.dem;
    
    ifstream ampfile(mlcfile.c_str(), ios::in | ios::binary);
    if (!ampfile.is_open()){
      cout << "Error opening input intensity file " << mlcfile << "\n";
        exit(1);
          }
    else
      cout << "Opened input intensity file: " << mlcfile << endl;
    
    ifstream DEMfile(hgtfile.c_str(), ios::in | ios::binary);
    if (!DEMfile.is_open()){
      cout << "Error opening DEM file " << hgtfile << "\n";
        exit(1);
          }
    else
      cout << "Opened DEM file: " << hgtfile << endl;
    
    if (area_flag){
      areaRDCout.open(area_out.c_str(), ios::out | ios::binary);
      if (!areaRDCout.is_open()){
        cout << "Error creating output area file " << area_out << "\n";
          exit(1);
        }
      else
        cout << "Created output area file: " << area_out << endl;
    }
    if (LUTin_flag){
      LUTin.open(LUT_flin.c_str(), ios::in | ios::binary);
      if (!LUTin.is_open()){
        cout << "Error opening input look up table " << LUT_flin << "\n";
          exit(1);
        }
      else
        cout << "Opened input look up table: " << LUT_flin << endl;
    }
    if (LUTout_flag){
      LUTout.open(LUT_flout.c_str(), ios::out | ios::binary);
      if (!LUTout.is_open()){
        cout << "Error creating output look up table " << LUT_flout << "\n";
          exit(1);
        }
      else
        cout << "Created output look up table: " << LUT_flout << endl;
    }
    if (sim_flag){
      sim_flout.open(sim_name.c_str(), ios::out | ios::binary);
      if (!sim_flout.is_open()){
        cout << "Error creating output simulated SAR image " << sim_name << "\n";
          exit(1);
        }
      else
        cout << "Created output simulated SAR image: " << sim_name << endl;
    }
    if (correct_flag){
      ampout.open(par.cor_out.c_str(), ios::out | ios::binary);
      if (!ampout.is_open()){
        cout << "Error creating output corrected intensity file " << par.cor_out << "\n";
          exit(1);
        }
      else
        cout << "Created output corrected intensity file: " << par.cor_out << endl;
    }
    if (cos_flag) {
      VegTablefile.open(veg_in.c_str(), ios::in | ios::binary);
      if (!VegTablefile.is_open()){
        cout << "Error opening input Vegetation Correction file " << veg_in << "\n";
        exit(1);
      }
      else
        cout << "Opened input Vegetation Correction file: " << veg_in << endl;
    }
    if (look_flag){
        look_out.open(look_name.c_str(), ios::out | ios::binary);
        if (!look_out.is_open()){
            cout << "Error creating output look angle file " << look_name << "\n";
            exit(1);
        }
        else
            cout << "Created output look angle file: " << look_name << endl;
    }
    if (slope_flag){
        slope_out.open(slope_name.c_str(), ios::out | ios::binary);
        if (!slope_out.is_open()){
            cout << "Error creating output slope angle file " << slope_name << "\n";
            exit(1);
        }
        else
            cout << "Created output slope angle file: " << slope_name << endl;
    }
    if (mask_flag){
        mask_out.open(mask_name.c_str(), ios::out | ios::binary);
        if (!mask_out.is_open()){
            cout << "Error creating output validity mask file " << mask_name << "\n";
            exit(1);
        }
        else
            cout << "Created output validity mask file: " << mask_name << endl;
    }
    if (ratio_flag){
        ratio_out.open(diff_name.c_str(), ios::out | ios::binary);
        if (!ratio_out.is_open()){
            cout << "Error creating output calibration difference file " << diff_name << "\n";
            exit(1);
        }
        else
            cout << "Created output calibration difference file: " << diff_name << endl;
    }

        
    //Create buffer vectors for data
    vector<float> zero_vec(par.widthDEM,0), simsar(par.widthDEM,0), gc1(par.widthDEM,0), gc2(par.widthDEM,0);
    vector<float> look_array(par.widthDEM,0), slope_array(par.widthDEM,0);
    vector<complex<float> > gc(par.widthDEM,0), gc_out(par.widthDEM,0), zero_vec_cpx(par.widthDEM,0);
    vector<vector<float> > DEM_buf_float(3,vector<float>(par.widthDEM,0));
    vector<vector<float> > VegTable( 900, vector<float> (900,0.0001) );

    //Create buffer vectors for JPL areas in RDC coordinates
    vector<float> area_fe_vec(par.width,0), diff_area_fe_vec(par.width,0);
        
    //Create 2-D vector-vectors for area and local incidence angle estimates in RDC coordinates
    vector<vector<float> > areaRDC(par.height,vector<float>(par.width,0)), theta_l(par.height,vector<float>(par.width,0)), 
      count(par.height,vector<float>(par.width,0)), sum_wgt(par.height,vector<float>(par.width,0)),
      r_looks(par.height,vector<float>(par.width,0)),all_slope_actual_r(par.height,vector<float>(par.width,0)), antcors(par.height,vector<float>(par.width,0));


    // --------------------   Main code: decompose DEM into facets, compute RDC coordinates, and area/local_inc in RDC  --------------------------

    //Estimate parameters at peg point
    peg.pos = llh2ecef(peg.lat, peg.lon, 0.0f);
    peg.re = WGS84_A/sqrt(1.0-WGS84_E2*sin(peg.lat)*sin(peg.lat));
    peg.rn = WGS84_A*(1.0-WGS84_E2)/sqrt(POW3(1.0-WGS84_E2*sin(peg.lat)*sin(peg.lat)));
    peg.ra = peg.re*peg.rn/(peg.re*cos(peg.heading)*cos(peg.heading) + peg.rn*sin(peg.heading)*sin(peg.heading)); //Radius of approximating sphere at peg
    temp_raU.x = 0; temp_raU.y = 0; temp_raU.z = peg.ra;
    peg.raU = ECEF_transform(temp_raU, peg.lat, peg.lon);
    r_x1 = peg.ra + par.gavgalt;
    r_x2 = peg.ra + par.gavgterhgt;
    
    //Estimate parameters and DEM spacing (in SCH sense) at image center
    lat = par.corner_lat - 0.5*par.heightDEM*par.spc_lat;
    lon = par.corner_lon + 0.5*par.widthDEM*par.spc_lon;
    SCH = llh2sch(lat, lon, 0.0, par, peg);
    lon = lon + par.spc_lon;
    SCH2 = llh2sch(lat, lon, 0.0, par, peg);
    deltaDEM_lon = sqrt(POW2(SCH2.x-SCH.x) + POW2(SCH2.y-SCH.y));
    lon = par.corner_lon + 0.5*par.widthDEM*par.spc_lon;
    lat = lat - par.spc_lat;
    SCH2 = llh2sch(lat, lon, 0.0, par, peg);
    deltaDEM_lat = sqrt(POW2(SCH2.x-SCH.x) + POW2(SCH2.y-SCH.y));
    area_ref = par.delta_az*par.delta_R; //reference RDC area per pixel



    //Enter loop to read or compute transformation LUT and area
    cout << "\n";
    i_cur = 0.0;
    for (long ii = 0; ii < par.heightDEM; ++ii){

        if (ii % 1000 == 0)
            cout << "Processed line " << ii << " of " << par.heightDEM << "\r" << flush;

        //If input LUT is provided, read values (cpx format/BIP)
        if (LUTin_flag)
            LUTin.read((char *) &gc[0], sizeof(float)*2*par.widthDEM);
                
        //Load 3-line input DEM buffer
        if (ii == 0 || ii == (par.heightDEM-1)){ //Check bounds
            if (sim_flag)
                sim_flout.write((char *) &zero_vec[0], sizeof(float)*par.widthDEM);
            if (LUTout_flag)
                LUTout.write((char *) &zero_vec_cpx[0], sizeof(float)*2*(par.widthDEM));
            if (look_flag)
                look_out.write((char *) &zero_vec[0], sizeof(float)*par.widthDEM);
            if (slope_flag)
                slope_out.write((char *) &zero_vec[0], sizeof(float)*par.widthDEM);
            // if (ratio_flag)
            //  ratio_out.write((char *) &zero_vec[0], sizeof(float)*par.widthDEM);

            i_cur += 1.0;
            continue;
        }
        else {
            DEMfile.seekg(sizeof(float)*(par.widthDEM*(ii-1)), ios::beg);
            for (short i = 0; i < 3; ++i)
                    DEMfile.read((char *) &DEM_buf_float[i][0], sizeof(float)*par.widthDEM);
        }
                
        lat = par.corner_lat - i_cur*par.spc_lat;
        j_cur = 0.0;


        for (long jj = 0; jj < par.widthDEM; ++jj){

            if (jj == 0 || jj == (par.widthDEM-1)){ //Check bounds
                simsar[jj] = 0.0f;
                j_cur += 1.0;
                continue;
            }

            //Load DEM into 3x3 window
            Z1 = DEM_buf_float[0][jj-1];
            Z2 = DEM_buf_float[0][jj];
            Z3 = DEM_buf_float[0][jj+1];
            Z4 = DEM_buf_float[1][jj-1];
            Z5 = DEM_buf_float[1][jj];
            Z6 = DEM_buf_float[1][jj+1];
            Z7 = DEM_buf_float[2][jj-1];
            Z8 = DEM_buf_float[2][jj];
            Z9 = DEM_buf_float[2][jj+1];

            //Z1 = par.gavgterhgt; Z2 = Z1; Z3 = Z1; Z4 = Z1; Z5 = Z1; Z6 = Z1; Z7 = Z1; Z8 = Z1; Z9 = Z1;

            if (LUTin_flag){
                ranpix = gc[jj].real();
                azpix = gc[jj].imag();
                slt_range = par.Ro + (double)ranpix*par.delta_R;
            }
            else {
                if (Z5 < -1000){ //Discard bad DEM data point (mainly for UAVSAR DEMs)
                    azpix = -100.0f;
                    ranpix = -100.0f;
                }
                else {
                    lon = par.corner_lon + j_cur*par.spc_lon;
                                        
                    //Compute SCH coordinates for map pixel
                    SCH = llh2sch(lat, lon, Z5, par, peg);
                    
                    //Convert SCH coordinates to range and azimuth positions
                    if (SCH.y < 0.0f){
                        azpix = -100.0f;
                        ranpix = -100.0f;
                    }
                    else {
                        azpix = (SCH.x-par.so)/par.delta_az;
                        slt_range = sqrt( POW2(peg.ra+(double)Z5) + POW2(peg.ra+par.gavgalt) - 2.0*(peg.ra+(double)Z5) * (peg.ra+par.gavgalt) * cos(SCH.y/peg.ra));
                        ranpix = (slt_range-par.Ro)/par.delta_R;
                    }
                }               
                if (LUTout_flag){
                  //gc_out[jj].real() = ranpix;
                  //gc_out[jj].imag() = azpix;
                  gc_out[jj] = std::complex<float>(ranpix,azpix);
     
                }
            }

            //Establish bounds for bilinear weighting model
            x1 = floor(ranpix); ix1 = (int)x1;
            x2 = x1+1.0f; ix2 = (int)x2;
            y1 = floor(azpix); iy1 = (int)y1;
            y2 = y1+1.0f; iy2 = (int)y2;        
                
            //Check to see if pixel lies in valid RDC range
            if (ranpix < 0.0f || x2 > xbound || azpix < 0.0f || y2 > ybound ){
                simsar[jj] = 0.0f;
                j_cur += 1.0;
                continue;
            }

            //Compute slope and aspect for pixel using 3x3 window
            p = (Z3 + Z6 + Z9 - Z1 - Z4 - Z7) / (6.0f * deltaDEM_lon);
                q = (Z1 + Z2 + Z3 - Z7 - Z8 - Z9) / (6.0f * deltaDEM_lat);
            slope = atan(sqrt(p*p + q*q));
            if (p == 0.0f){
                if (q > 0)
                    aspect = PI;
                else
                    aspect = 0.0f;
            }
            else
                aspect = PI - atan(q/p) + (PI_HALF*p/fabs(p));

            //Slope in the range and azimuth directions (for left-looking sensor)
            slope_r = tan(slope)*cos(aspect - peg.heading - PI_HALF); // (-PI_HALF) indicates slopes towards radar are positive
            slope_a = tan(slope)*cos(aspect - peg.heading);  //THESE values are tan(actual_slope_a)
            slope_actual_r = atan(slope_r);
            slope_actual_a = atan(slope_a);
            //Vector (sch) of unit normal vector to surface adjusted for range pixel
            temp = -1.0/sqrt(1.0 + POW2(slope_r) + POW2(slope_a));
            nE.x = temp*slope_a;
            nE.y = temp*slope_r;
            nE.z = -temp;           

            //Compute look vector (direction from ground to sensor ---> negative of the convention)
            r_x2 = peg.ra+Z5;
            r_l3 = (r_x1*r_x1 + slt_range*slt_range - r_x2*r_x2)/(2.0*r_x1*slt_range);
            r_look = acos((r_l3 + sin(par.ESA)*sin(par.pitch))/(cos(par.pitch)*cos(par.ESA)));
            
            // added by Michael Denbina to put look and slope into arrays for saving.
            look_array[jj] = r_look*(180.0/PI);
            slope_array[jj] = slope_actual_r*(180.0/PI);


            nL.x = -sin(par.ESA)*cos(par.pitch)*cos(par.yaw)
                   -cos(par.ESA)*(sin(par.pitch)*cos(r_look)*cos(par.yaw) + sin(r_look)*sin(par.yaw));
            nL.y = sin(par.ESA)*cos(par.pitch)*sin(par.yaw) 
                  -cos(par.ESA)*(-sin(par.pitch)*cos(r_look)*sin(par.yaw) + sin(r_look)*cos(par.yaw));
            nL.z = r_l3;

            //Compute normal vector to imaging plane (remember the reverse direction)
            nI.x = 0.0;
            nI.y = nL.z;
            nI.z = -nL.y;

            //Compute local incidence angle and map area for facet
            area = area_ref/fabs(dotXYZ(nE,nI));
            temp = acos(dotXYZ(nE,nL)); //local incidence angle
            if (sim_flag){
              simsar[jj] =  temp; //atan(slope_a); //temp;  r_look;
            }

            // two-way amplitude gain of antenna
            //antcor = pow(sin((acos(par.gavgalt/slt_range)-0.785398)*PI*1.5)/((acos(par.gavgalt/slt_range)-0.785398)*PI*1.5),2);
            //antcor = pow(sin((acos((par.gavgalt-par.gavgterhgt)/slt_range)-0.785398)*PI*1.5)/((acos((par.gavgalt-par.gavgterhgt)/slt_range)-0.785398)*PI*1.5),2);
            // remove the correction by using actual look angle instead of approximation makes antcor equal to 1
            antcor = pow(sin((r_look-0.785398)*PI*1.5)/((r_look-0.785398)*PI*1.5),2);

            //Use 1/dist IDW to assign incidence angle and area to radar pixel
            dist = sqrt(POW2(x1-ranpix) + POW2(y1-azpix));
            areaRDC[iy1][ix1] += area/dist;
            theta_l[iy1][ix1] += temp/dist;
            sum_wgt[iy1][ix1] += 1.0f/dist;
            r_looks[iy1][ix1] += r_look/dist;
            all_slope_actual_r[iy1][ix1] += slope_actual_r/dist;
            antcors[iy1][ix1] += antcor/dist;

            dist = sqrt(POW2(x2-ranpix) + POW2(y1-azpix));
            areaRDC[iy1][ix2] += area/dist;
            theta_l[iy1][ix2] += temp/dist;
            sum_wgt[iy1][ix2] += 1.0f/dist;
            r_looks[iy1][ix2] += r_look/dist;
            all_slope_actual_r[iy1][ix2] += slope_actual_r/dist;
            antcors[iy1][ix2] += antcor/dist;

            dist = sqrt(POW2(x1-ranpix) + POW2(y2-azpix));
            areaRDC[iy2][ix1] += area/dist;
            theta_l[iy2][ix1] += temp/dist;
            sum_wgt[iy2][ix1] += 1.0f/dist;
            r_looks[iy2][ix1] += r_look/dist;
            all_slope_actual_r[iy2][ix1] += slope_actual_r/dist;
            antcors[iy2][ix1] += antcor/dist;

            dist = sqrt(POW2(x2-ranpix) + POW2(y2-azpix));
            areaRDC[iy2][ix2] += area/dist;
            theta_l[iy2][ix2] += temp/dist;
            sum_wgt[iy2][ix2] += 1.0f/dist;
            r_looks[iy2][ix2] += r_look/dist;
            all_slope_actual_r[iy2][ix2] += slope_actual_r/dist;
            antcors[iy2][ix2] += antcor/dist;
            
            //r_looks[ii][jj] = r_look;
            //all_slope_actual_r[ii][jj] = slope_actual_r;

            j_cur += 1.0; //Update counter

        }
        if (LUTout_flag)
            LUTout.write((char *) &gc_out[0], sizeof(float)*2*(par.widthDEM));
        
        if (look_flag)
            look_out.write((char *) &look_array[0], sizeof(float)*par.widthDEM);
        
        if (slope_flag)
            slope_out.write((char *) &slope_array[0], sizeof(float)*par.widthDEM);
        
        if (sim_flag)
            sim_flout.write((char *) &simsar[0], sizeof(float)*(par.widthDEM));
            
        i_cur += 1.0; //Update counter
    }



    //Compute weighted incidence angle and area matrices
    for (long i = 0; i < par.height; ++i){
        for (long j = 0; j < par.width; ++j){
            if (sum_wgt[i][j] < 1.0e-4){
                theta_l[i][j] = 0.0f;
                areaRDC[i][j] = 0.0f;
                r_looks[i][j] = 0.0f;
                all_slope_actual_r[i][j] = 0.0f;
                antcors[i][j] =0.0f;

                continue;
            }
            theta_l[i][j] /= sum_wgt[i][j];
            areaRDC[i][j] /= sum_wgt[i][j];
            r_looks[i][j] /= sum_wgt[i][j];
            all_slope_actual_r[i][j] /= sum_wgt[i][j];
            antcors[i][j] /= sum_wgt[i][j];

        }
    }
    
    //(Optional) Write out RDC area estimate to file
    if (area_flag){
        for (long i = 0; i < par.height; ++i)
            areaRDCout.write((char *) &areaRDC[i][0], sizeof(float)*par.width);
        areaRDCout.close();
    }

    /*//Write out incidence angle in RDC coordinates
    inc_tol = 87.0*RAD;
    cout << "Writing out incidence angle image........" << flush;
    ofstream inc_flout("slope_range_rdc.bin", ios::out | ios::binary);
    for (long i = 0; i < par.height; ++i){
        inc_flout.write((char *) &theta_l[i][0], sizeof(float)*par.width);
    }
    inc_flout.close();
    cout << "Done\n\n";*/
    
    //-----------------------------   Optional: Perform correction to intensity image  -----------------------------------

    if (correct_flag){
      
      cout << "\n\nCorrecting input intensity image " << flush;
      
      vector<float> amp_in(par.width,0), amp_cor(par.width,0), mask_array(par.width,0), rtc_ratio(par.width,0);

        //Compute 1-D look-up vectors for JPL area correction factors
        compute_area_fe(peg, par, area_fe_vec);
        
        //Enter loop to read data and perform radiometric correction
        inc_tol = 70.0*RAD; //tolerance for local incidence angle in shadows
        inc_ltol= 15*RAD;  // low incidence angle. data gets strectched
        if (cos_flag){
            cout << "using Vegetation correction ....." << flush;
            // reading and processing backscatter data
              cout << "Reading Vegetation Correction table ....." << flush;
              for (short i = 0; i < 900; ++i)
                VegTablefile.read((char *) &VegTable[i][0], sizeof(float)*900);

            for (long i = 0; i < par.height; ++i){
                ampfile.read((char *) &amp_in[0], sizeof(float)*par.width);         
                for (long j = 0; j < par.width; ++j) {
                    // Preserve original values
                    //amp_og[j] = amp_in[j];
                  
                    // Remove JPL correction factor
                    //amp_in[j] = amp_in[j]*area_fe_vec[j];
                    rtc_ratio[j] = area_fe_vec[j];

                    if (areaRDC[i][j] < 1.0e-6) { // Negligible area calculated for coordinate
                      //amp_cor[j] = (amp_in[j]);
                      mask_array[j] = 1; // Added by Michael Denbina to keep track of void pixels. 1 = void pixel
                    }
                    else {
                        inc_cor = theta_l[i][j];
                        r_look = r_looks[i][j];
                        slope_actual_r=all_slope_actual_r[i][j];

                        //tempout = amp_in[j]*(area_ref/areaRDC[i][j])/cos(inc_cor);  // AREA CORRECTION
                        antcor = antcors[i][j] / 
                          pow(sin((r_look-0.785398)*PI*1.5)/
                              ((r_look-0.785398)*PI*1.5),2); // antenna correction due to original image not using DEM
                        //tempout = tempout * antcor;  // ANTENNA SUPPLEMENTAL CORRECTION FOR TOPO
                        rtc_ratio[j] = rtc_ratio[j]*(area_ref/areaRDC[i][j])/(cos(inc_cor))*antcor;

                        if ( (inc_cor > inc_tol) || (inc_cor < inc_ltol) || (r_look < 0.35 && abs(slope_actual_r) > 0.1) ) {
                          //amp_cor[j] = tempout;
                          mask_array[j] = 1; // Added by Michael Denbina to keep track of void pixels. 1 = void pixel
                        }
                        else {
                          mask_array[j] = 0; // Added by Michael Denbina to keep track of void pixels.  0 = valid pixel


                          e_look = (int) (r_look*180/PI)*10; 
                          e_slope= (int) (slope_actual_r*180/PI+90.)*10/2;
                          
                          if(e_look < 0 || e_look > 899 || e_slope <0 || e_slope >899 ) cs = 1.0;
                          else cs = VegTable[e_slope][e_look]; // VEGETATION Table lookup
                          
                          // Testing:


                          if (cs > 0.001)
                          switch (pol) {
                          case 1: 

                            //cs = 0.816471 -0.933689*r_look +0.251198 *pow(r_look,2)+0.246271*slope_actual_r-0.236941*pow(slope_actual_r,2);
                            //cs = 2.38234 -2.69612*r_look+ 0.709170*pow(r_look,2)+0.726989*slope_actual_r -0.700571*pow(slope_actual_r,2);
                            
                            // cout<<" look "<<e_look<<" slope "<<e_slope<<" ";
                            //if(cs = 0) amp_cor[j]=1;
                            //amp_cor[j] = (tempout/cs*VegTable[450][350]);
                            rtc_ratio[j] = rtc_ratio[j]/cs*VegTable[450][350];
                            break;
                            
                          case 2: 
                            
                            //cs = 0.255000-0.311589 *r_look +0.0980392*pow(r_look,2)+0.0621085*slope_actual_r -0.0633426*pow(slope_actual_r,2);
                            //cs = 2.49554 -3.01881*r_look+ 0.934043 *pow(r_look,2)+ 0.615928*slope_actual_r -0.629401*pow(slope_actual_r,2);
                            
                            //amp_cor[j] = (tempout/cs*VegTable[450][350]);
                            rtc_ratio[j] = rtc_ratio[j]/cs*VegTable[450][350];
                            break;
                            
                          case 3:

                            //cs= 0.673232 -0.922350 *r_look +0.343811*pow(r_look,2)+0.128466*slope_actual_r -0.125623*pow(slope_actual_r,2);
                            //cs= 2.79799 -3.79865*r_look+1.40014*pow(r_look,2)+0.542621*slope_actual_r-0.532302*pow(slope_actual_r,2);
                            
                            //amp_cor[j] = (tempout/cs*VegTable[450][350]);
                            rtc_ratio[j] = rtc_ratio[j]/cs*VegTable[450][350];
                            break;
                            
                            
                          // default:
                            
                          //   amp_cor[j] = tempout/cs;
                          }
                          else {
                              mask_array[j] = 1;
                              //amp_cor[j] = tempout;
                          }
                        }
                    }
                    if (!(amp_cor[j] <= DBL_MAX && amp_cor[j] >= -DBL_MAX)) {
                        //amp_cor[j] = 0;
                        mask_array[j] = 1;
                    }
                    
                    // Difference
                    //rtc_ratio[j] = amp_cor[j]/amp_og[j];
                    
                    // Apply upper and lower limits to correction factor.
                    if (rtc_ratio[j] < min_correction_ratio) {
                        rtc_ratio[j] = min_correction_ratio;
                    }
                    else if (rtc_ratio[j] > max_correction_ratio) {
                        rtc_ratio[j] = max_correction_ratio;
                    }
                    
                    // If mask is set to void, set correction factor to void.
                    // and keep input value as is
                    if (mask_array[j] == 1) {
                        rtc_ratio[j] = void_correction_val;
                        //amp_cor[j] = amp_in[j];
                        amp_cor[j] = void_correction_val;
                    }
                    
                    // Update amp_cor, in case rtc_ratio changed.
                    if (rtc_ratio[j] > 0) {
                        amp_cor[j] = rtc_ratio[j] * amp_in[j];
                    } 
                }           
                //Write out corrected data in RDC coordinates
                ampout.write((char *) &amp_cor[0], sizeof(float)*par.width);

                if (ratio_flag)
                    ratio_out.write((char *) &rtc_ratio[0], sizeof(float)*par.width);
                
                if (mask_flag)
                    mask_out.write((char *) &mask_array[0], sizeof(float)*par.width);
            }
        }
        else {
            cout << "using area correction....." << flush;
            for (long i = 0; i < par.height; ++i){
                ampfile.read((char *) &amp_in[0], sizeof(float)*par.width);         
                for (long j = 0; j < par.width; ++j) {
                    ////Preserve original values
                    //amp_og[j] = amp_in[j];
                    ////Remove JPL correction factor
                    //amp_in[j] = amp_in[j]*area_fe_vec[j];
                    r_look = r_looks[i][j];
                    antcor = antcors[i][j] / 
                             pow(sin((r_look-0.785398)*PI*1.5)/
                             ((r_look-0.785398)*PI*1.5),2); // antenna correction due to original image not using DEM
                    //amp_in[j] = amp_in[j]* antcor;
                    rtc_ratio[j] = area_fe_vec[j]*antcor;
                                        
                    if (areaRDC[i][j] < 1.0e-6) { //No area calculated for coordinate
                        //tempout = amp_in[j];
                        mask_array[j] = 1; // void pixel
                    }
                    else {
                        inc_cor = theta_l[i][j];
                        if (inc_cor > inc_tol || (inc_cor < inc_ltol)) { // || (r_look < 0.35 && abs(slope_actual_r) > 0.1)) //Pixel is in shadow
                          //tempout = amp_in[j];
                          mask_array[j] = 1; // void pixel
                        }
                        else {
                          //tempout = amp_in[j]*(area_ref/areaRDC[i][j])/cos(inc_cor);
                          rtc_ratio[j] = rtc_ratio[j]*(area_ref/areaRDC[i][j])/cos(inc_cor);
                          mask_array[j] = 0; // valid pixel
                        }

                    }
                    //amp_cor[j] = tempout;
                    //rtc_ratio[j] = amp_cor[j]/amp_og[j];
                    amp_cor[j] = rtc_ratio[j] * amp_in[j];

                    // Apply upper and lower limits to correction factor.
                    if (rtc_ratio[j] < min_correction_ratio) {
                        rtc_ratio[j] = min_correction_ratio;
                    }
                    else if (rtc_ratio[j] > max_correction_ratio) {
                        rtc_ratio[j] = max_correction_ratio;
                    }
                    
                    // If mask is set to void, set correction factor to void.
                    // and keep input value as is
                    if (mask_array[j] == 1) {
                        rtc_ratio[j] = void_correction_val;
                        //amp_cor[j] = amp_in[j];
                        amp_cor[j] = void_correction_val;
                    }
                    
                    // Update amp_cor, in case rtc_ratio changed.
                    if (rtc_ratio[j] > 0) {
                        amp_cor[j] = rtc_ratio[j] * amp_in[j];
                    } 
                }           
                //Write out corrected data in RDC coordinates
                ampout.write((char *) &amp_cor[0], sizeof(float)*par.width);
                
                if (mask_flag)
                    mask_out.write((char *) &mask_array[0], sizeof(float)*par.width);

                if (ratio_flag)
                    ratio_out.write((char *) &rtc_ratio[0], sizeof(float)*par.width);
            }
        }
        cout << "Done" << endl;
        amp_cor.clear(); amp_in.clear(); rtc_ratio.clear(); 
    }


    DEM_buf_float.clear(); theta_l.clear(); areaRDC.clear(); zero_vec.clear(); LUTcpx.clear(); count.clear(); gc.clear(); gc_out.clear(); simsar.clear();r_looks.clear();all_slope_actual_r.clear(); look_array.clear(); slope_array.clear();
    DEMfile.close(); LUTout.close(); ampfile.close(); ampout.close(); areaRDCout.close(); sim_flout.close();
    
    if (ratio_flag)
        ratio_out.close();

    if (mask_flag)
        mask_out.close();
    
    if (look_flag)
        look_out.close();
    
    if (slope_flag)
        slope_out.close();
    
    timerend = time (NULL); //End timer
    cout << "\nElapsed time: " << timerend-timerstart << " seconds" << endl;

    return 0;
}






    
