#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <string>
#include <math.h>
#include <cstring>
#include <cstdlib>
using namespace std;

int main(int argc, char* argv[]){

	int ix1, ix2, iy1, iy2, i_bound_first_flag;
	long width, height, size, i_bound_first, i_bound_last, j_bound_first, j_bound_last, width_LUT, height_LUT, i_stop, j_stop, j_out;
	float x1, x2, y1, y2, denom, ranpix, azpix, xbound, ybound;
	double corner_lat, corner_lon, temp, deg_unit;
	string name_int, name_LUT, name_out, lutrsc_name, rsc_name;

	//Get command line arguments
	switch (argc){
		case 7:
			name_int = argv[1];
			width = atol(argv[2]);
			name_LUT = argv[3];
			name_out = argv[4];
			width_LUT = atol(argv[5]);
			height_LUT = atol(argv[6]);
			break;
		default:
			cout << "Error: invalid number of arguments\n";
			cout << "Usage: " << argv[0] << " <name_int>  <width in>  <LUT>  <name_out>  <width out>  <height out> \n";
			exit(1);
	}

	cout << "Geocoding intensity files from RDC to MAP coordinates\n";
	cout << "Input file: " << name_int << "\n";
	cout << "Output file: " << name_out << "\n";

	//Create file pointers
	ifstream int_flin(name_int.c_str(), ios::in | ios::binary);
	if (!int_flin.is_open()){
		cout << "Error opening intensity image\n";
		exit(1);
	}
	ifstream LUT_flin(name_LUT.c_str(), ios::in | ios::binary);
	if (!LUT_flin.is_open()){
		cout << "Error opening transformation LUT\n";
		exit(1);
	}
	ofstream int_flout(name_out.c_str(), ios::out | ios::binary);
	if (!int_flout.is_open()){
		cout << "Error creating geocoded intensity image\n";
		exit(1);
	}

	//Determine number of lines in intensity image
	int_flin.seekg(0, ios::end);
	size = int_flin.tellg();
	height = size/(sizeof(float)*width);
	int_flin.seekg(0, ios::beg);
	xbound = (float)width-1.0f;
	ybound = (float)height-1.0f;

	//Create vector-vectors for data
	vector<vector<float> > int_in(height, vector<float>(width,0));
	vector<complex<float> > LUT_in(width_LUT,0);
	vector<float> int_out(width_LUT,0);
	
	//Load input intensity data
	for (long i = 0; i < height; ++i)
		int_flin.read((char *) &int_in[i][0], sizeof(float)*width);

	cout << "\n";
	//Enter loop to geocode data using bilinear interpolation
	for (long i = 0; i < height_LUT; ++i){

		if (i % 1000 == 0)
			cout << "Processed line " << i << " of " << height_LUT << "\r" << flush;
		
		//Read in CPX format LUT
		LUT_flin.read((char *) &LUT_in[0], sizeof(float)*2*width_LUT);
		
		for (long j = 0; j < width_LUT; ++j){

			ranpix = LUT_in[j].real();
			azpix = LUT_in[j].imag();

			if (ranpix <= 0 || ranpix >= xbound || azpix <= 0 || azpix >= ybound){
				int_out[j] = 0.0f;
				continue;
			}

			x1 = floor(ranpix); ix1 = (int)x1;
			x2 = ceil(ranpix); ix2 = (int)x2;
			y1 = floor(azpix); iy1 = (int)y1;
			y2 = ceil(azpix); iy2 = (int)y2;
			if (fabs(ranpix-x1) < 1.0e-5 && fabs(azpix-y1) < 1.0e-5)
				int_out[j] = int_in[iy1][ix1];
			else if (fabs(ranpix-x1) < 1.0e-5)
				int_out[j] = (y2-azpix)*int_in[iy1][ix1] + (azpix-y1)*int_in[iy2][ix1];
			else if (fabs(azpix-y1) < 1.0e-5)
				int_out[j] = (x2-ranpix)*int_in[iy1][ix1] + (ranpix-x1)*int_in[iy1][ix2];
			else {						
				int_out[j] = (x2-ranpix)*(y2-azpix)*int_in[iy1][ix1]
				       + (ranpix-x1)*(y2-azpix)*int_in[iy1][ix2]
   				       + (x2-ranpix)*(azpix-y1)*int_in[iy2][ix1]
				       + (ranpix-x1)*(azpix-y1)*int_in[iy2][ix2];
			}
						
		}
		
		int_flout.write((char *) &int_out[0], sizeof(float)*width_LUT);
	}

	cout << "\n\nDone" << endl;

	int_flin.close(); LUT_flin.close(); int_flout.close();
	int_in.clear(); LUT_in.clear(); int_out.clear();

	return 0;

}

