#ifndef FACET_MATH_UAVSAR_H
#define FACET_MATH_UAVSAR_H

#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>
#include <vector>

#define PI     3.141592653589793
#define PI_HALF     1.5707963267949
#define RAD     0.0174532925199433
#define WGS84_A     6378137.0
#define WGS84_E2    0.00669437999015

#define POW2(a)     ((a)*(a))
#define POW3(a)     ((a)*(a)*(a))
#define POW4(a)     ((a)*(a)*(a)*(a))
#define POW5(a)     ((a)*(a)*(a)*(a)*(a))
#define POW6(a)     ((a)*(a)*(a)*(a)*(a)*(a))
#define POW7(a)     ((a)*(a)*(a)*(a)*(a)*(a)*(a))
#define POW8(a)     ((a)*(a)*(a)*(a)*(a)*(a)*(a)*(a))

struct par_struct {	
	int rlks, azlks;
	long width, height, widthDEM, heightDEM;
	double delta_az, delta_R, delta_t_az, Ro, so, co, spc_lat, spc_lon, corner_lat, corner_lon, glob_inc, az, gavgalt, gavgterhgt, 
               pitch, ESA, yaw, d;
	std::string ann, mlc, dem, cor_out, pol;
};

struct XYZ {
	double x, y, z;
};

struct peg_struct{

	double lat, lon, heading, re, rn, ra;
	XYZ pos, raU;

};

struct LLH_info {
	double lat, lon, h;
};

XYZ addXYZ (XYZ a, XYZ b){

	XYZ c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return c;
}

XYZ subXYZ (XYZ a, XYZ b){

	XYZ c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return c;
}

XYZ crossXYZ(XYZ a, XYZ b){

	XYZ c;
	c.x = a.y*b.z - a.z*b.y;
	c.y = -(a.x*b.z - a.z*b.x);
	c.z = a.x*b.y - a.y*b.x;
	return c;
}

double normXYZ (XYZ a){

	double result;
	result = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
	return result;
}

double dotXYZ (XYZ a, XYZ b){

	double product;
	product = a.x*b.x + a.y*b.y + a.z*b.z;
	return product;
}

int round_number(double param){

	double intpart, fractpart, rounded;
	int round_int;

	fractpart = modf(param, &intpart);
	if (fractpart >= 0.5) {
		rounded = ceil(param);
	}
	else {
		rounded = floor(param);
	}
	round_int = int(rounded);
	return round_int;
}

float fround_number(double param){

	double intpart, fractpart, rounded;
	float round_float;

	fractpart = modf(param, &intpart);
	if (fractpart >= 0.5) {
		rounded = ceil(param);
	}
	else {
		rounded = floor(param);
	}
	round_float = float(rounded);
	return round_float;
}

XYZ llh2ecef(double lat, double lon, float h){

	//This subroutine converts latitude/longitude/height to ECEF XYZ coordinates
	//Default values are for WGS-84 ellipsoid
	
	double Nh;
	XYZ X;

	//Ellipsoid parameters
	Nh = WGS84_A/sqrt(1.0-WGS84_E2*sin(lat)*sin(lat));

	//Convert to XYZ
	X.x = (Nh+(double)h)*cos(lat)*cos(lon);
	X.y = (Nh+(double)h)*cos(lat)*sin(lon);
	X.z = (Nh+(double)h-WGS84_E2*Nh)*sin(lat);

	return X;

}

XYZ ECEF_transform (XYZ rhoT, double lat, double lon){

	//This subroutine convert topocentric vectors to ECEF vectors

	XYZ rho;
	double clat, clon, slat, slon;

	clat = cos(lat);
	clon = cos(lon);
	slat = sin(lat);
	slon = sin(lon);

	rho.x = -slon*rhoT.x - slat*clon*rhoT.y + clat*clon*rhoT.z;
	rho.y = clon*rhoT.x - slat*slon*rhoT.y + clat*slon*rhoT.z;
	rho.z = clat*rhoT.y + slat*rhoT.z;
	return rho;
}

XYZ llh2sch(double lat, double lon, double h, par_struct par, peg_struct peg){

	//This subroutine converts latitude/longitude/height to SCH
	//Default values are for WGS-84 ellipsoid
	
	double Nh, Nh_peg, rn, ra, ceta, seta, clat, slat, clon, slon, clat_peg, slat_peg, clon_peg, slon_peg, clambda, stheta, 
               T11, T12, T13, T21, T22, T23, T31, T32, T33;

	XYZ X, UTC, SCH, temp;

	//Compute trigonometric operations for later use
	ceta = cos(peg.heading);
	seta = sin(peg.heading);
	clat = cos(lat);
	slat = sin(lat);
	clon = cos(lon);
	slon = sin(lon);
	clat_peg = cos(peg.lat);
	slat_peg = sin(peg.lat);
	clon_peg = cos(peg.lon);
	slon_peg = sin(peg.lon);
	
	//Ellipsoid parameters
	Nh = WGS84_A/sqrt(1.0-WGS84_E2*slat*slat);
	
	//Convert current DEM point to WGS84 XYZ
	X.x = (Nh+h)*clat*clon;
	X.y = (Nh+h)*clat*slon;
	X.z = (Nh+h-WGS84_E2*Nh)*slat;

	//std::cout << std::setw(20) << std::setprecision(15) << X.x << std::setw(20) << std::setprecision(15) << X.y << std::setw(20) << std::setprecision(15) << X.z << "\n";
	
	//Temp position vector
	temp = subXYZ(X, subXYZ(peg.pos, peg.raU));
		
	//Rotation matrix from ECEF to UTC
	T11 = clat_peg*clon_peg;
	T12 = clat_peg*slon_peg;
	T13 = slat_peg;
	T21 = -seta*slon_peg-ceta*clon_peg*slat_peg;
	T22 = clon_peg*seta-ceta*slat_peg*slon_peg;
	T23 = ceta*clat_peg;
	T31 = ceta*slon_peg-seta*clon_peg*slat_peg;
	T32 = -clon_peg*ceta-seta*slat_peg*slon_peg;
	T33 = clat_peg*seta;

	//Compute UTC coordinates
	UTC.x = T11*temp.x + T12*temp.y + T13*temp.z;
	UTC.y = T21*temp.x + T22*temp.y + T23*temp.z;
	UTC.z = T31*temp.x + T32*temp.y + T33*temp.z;
	
	//Compute SCH coordinates
	clambda = atan(UTC.z/sqrt(UTC.x*UTC.x + UTC.y*UTC.y));
	stheta = atan(UTC.y/UTC.x);
	SCH.x = peg.ra*stheta;
	SCH.y = peg.ra*clambda;
	SCH.z = h;

	return SCH;
	
}

void compute_area_fe(peg_struct peg, par_struct par, std::vector<float> &area){

	//This subroutine computes the area correction factor applied to UAVSAR images in order to remove it

	double slt_range, r_x1, r_x2, r_l3, r_look, r_sininc, area1, area2;

	for (long j = 0; j < par.width; ++j){

		//Slant range for current pixel
		slt_range = par.Ro + (double)j*par.delta_R;

		//Compute area assuming a non-zero ESA and pitch
		r_x1 = peg.ra + par.gavgalt;
		r_x2 = peg.ra + par.gavgterhgt;
		r_l3 = (r_x1*r_x1 + slt_range*slt_range - r_x2*r_x2)/(2.0*r_x1*slt_range);
		r_look = acos((r_l3 + sin(par.ESA)*sin(par.pitch))/(cos(par.pitch)*cos(par.ESA)));
		r_sininc = sin(r_look)*r_x1/r_x2;
		area[j] = 1.0f/r_sininc;  ///par.delta_az*par.delta_R/r_sininc;// depends if UAVSAR processor used a reference area
	}
}



#endif
