#ifndef LOAD_ANN_FACET_H
#define LOAD_ANN_FACET_H

#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "math_uavsar.h"

void load_ann(par_struct &par, peg_struct &peg){

	//This subroutine will load auxiliary information from UAVSAR annotation file *.ann
	
    std::string dum[6], look;
	double slc_width, R_ground, delta_R_slc;

    //Open annotation file and read necessary info
    std::string line;
	std::ifstream infofile;
	infofile.open(par.ann.c_str(), std::ifstream::in);
	if (!infofile.is_open()){
		std::cout << "Error opening annotation file\n";
		exit(1);
	}
	while (getline(infofile, line)){
		std::istringstream stream(line, std::istringstream::in);
		stream >> dum[0];
		if (dum[0] == "mlc" + par.pol){
			stream >> dum[1] >> par.mlc;
			continue;
		}
		else if (dum[0] == "mlc_pwr.set_rows"){
			stream >> dum[1] >> dum[2] >> par.height;
			continue;
		}
		else if (dum[0] == "mlc_pwr.set_cols"){
			stream >> dum[1] >> dum[2] >> par.width;
			continue;
		}
		else if (dum[0] == "mlc_pwr.row_addr"){
			stream >> dum[1] >> dum[2] >> par.so;
			continue;
		}
		else if (dum[0] == "mlc_pwr.col_addr"){
			stream >> dum[1] >> dum[2] >> par.co;
			continue;
		}
		else if (dum[0] == "slc_mag.col_mult"){
			stream >> dum[1] >> dum[2] >> delta_R_slc;
			continue;
		}
		else if (dum[0] == "mlc_pwr.row_mult"){
			stream >> dum[1] >> dum[2] >> par.delta_az;
			continue;
		}
		else if (dum[0] == "mlc_pwr.col_mult"){
			stream >> dum[1] >> dum[2] >> par.delta_R;
			continue;
		}
		else if (dum[0] == "hgt"){
			stream >> dum[1] >> par.dem;
			continue;
		}
		else if (dum[0] == "hgt.set_rows"){
			stream >> dum[1] >> dum[2] >> par.heightDEM;
			continue;
		}
		else if (dum[0] == "hgt.set_cols"){
			stream >> dum[1] >> dum[2] >> par.widthDEM;
			continue;
		}
		else if (dum[0] == "hgt.col_addr"){
			stream >> dum[1] >> dum[2] >> par.corner_lon;
			par.corner_lon = par.corner_lon*RAD;
			continue;
		}
		else if (dum[0] == "hgt.row_addr"){
			stream >> dum[1] >> dum[2] >> par.corner_lat;
			par.corner_lat = par.corner_lat*RAD;
			continue;
		}
		else if (dum[0] == "hgt.col_mult"){
			stream >> dum[1] >> dum[2] >> par.spc_lon;;
			par.spc_lon = fabs(par.spc_lon)*RAD;			
			continue;
		}
		else if (dum[0] == "hgt.row_mult"){
			stream >> dum[1] >> dum[2] >> par.spc_lat;;
			par.spc_lat = fabs(par.spc_lat)*RAD;			
			continue;
		}
		else if (dum[0] == "Peg"){
			stream >> dum[1];
			if (dum[1] == "Latitude"){
				stream >> dum[2] >> dum[3] >> peg.lat;
				peg.lat = peg.lat*RAD;
			}
			else if (dum[1] == "Longitude"){
				stream >> dum[2] >> dum[3] >> peg.lon;
				peg.lon = peg.lon*RAD;
			}
			else if (dum[1] == "Heading"){
				stream >> dum[2] >> dum[3] >> peg.heading;
				peg.heading = peg.heading*RAD;
			}
			continue;
		}
		else if (dum[0] == "Image"){
			stream >> dum[1] >> dum[2] >> dum[3] >> dum[4] >> par.Ro;
			par.Ro = par.Ro*1000;
			continue;
		}
		else if (dum[0] == "Look"){
			stream >> dum[1] >> dum[2] >> dum[3] >> look;
			continue;
		}
		else if (dum[0] == "Global"){
			stream >> dum[1] >> dum[2];
			if (dum[2] == "Yaw"){
				stream >> dum[3] >> dum[4] >> par.yaw;
				par.yaw = par.yaw*RAD;
				continue;
			}
			else if (dum[2] == "Pitch"){
				stream >> dum[3] >> dum[4] >> par.pitch;
				par.pitch = par.pitch*RAD;
				continue;
			}
			else if (dum[2] == "ESA"){
				stream >> dum[3] >> dum[4] >> par.ESA;
				par.ESA = par.ESA*RAD;
				continue;
			}
			else if (dum[2] == "Altitude"){
				stream >> dum[3] >> dum[4] >> par.gavgalt;
				continue;
			}
			else if (dum[2] == "Terrain"){
				stream >> dum[3] >> dum[4] >> dum[5] >> par.gavgterhgt;
				continue;
			}			
			else
				continue;
		}
	}
	par.Ro += 1.0*delta_R_slc; //new addition to help compensate for range shift; basically, moving initial range to middle of 3-look window

	//Compute azimuth based on look direction (still only a rough estimate)
	if (look == "Right"){
		std::cerr << "Error: look direction is not 'Left'; must change SCH code first; exiting\n";
		exit(1);
	}
	else if (look == "Left"){
		par.az = peg.heading + PI_HALF;
	}
	else {
		std::cout << "Error computing azimuth\n";
		exit(1);
	}
	
	infofile.close();
}

#endif


