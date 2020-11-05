/************************************************************************
 * $Id: sunpos.h 3351 2018-02-21 15:00:33Z tobias.koelling $
 ************************************************************************/

/* This file is available in electronic form at http://www.psa.es/sdg/sunpos.htm */

#ifndef __sunpos_h
#define __sunpos_h

/* Declaration of some constants  */
#define pi    3.14159265358979323846
#define twopi (2*pi)
#define rad   (pi/180)
#define dEarthMeanRadius     6371.01	/* In km */
#define dAstronomicalUnit    149597890	/* In km */

#ifdef __cplusplus
extern "C" {
#endif

struct cTime
{
	int iYear;
	int iMonth;
	int iDay;
	double dHours;
	double dMinutes;
	double dSeconds;
};

struct cLocation
{
	double dLongitude;
	double dLatitude;
};

struct cSunCoordinates
{
	double dZenithAngle;
	double dAzimuth;
};

void sunpos(struct cTime udtTime, struct  cLocation udtLocation, struct cSunCoordinates *udtSunCoordinates);

#ifdef __cplusplus
}
#endif

#endif

