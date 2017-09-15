#ifndef __PHYSCONST_H
#define __PHYSCONST_H

#include <string>
#include <math.h>
using namespace std;

class PhysConst {
    public : 
        // class identifiers
        string name; 
        string linestyle;
        string markerstyle;
        string colorstyle;
        string savefilename;
        // mathematical constants //
        double pi;
        double piby2; 
        double sqrt2;
        double ln2;
        // astronomical constants //
        double earthradius;
        double sunradius;
        ///// physics constants/////
        double GF;
        double Na;
        double sw_sq;
        double G;
        double alpha;
        /////////// units //////////
        // energy
        double TeV;
        double GeV;
        double MeV;
        double keV;
        double Joule;
        double erg;
        // mass
        double kg;
        double gr;
        // time
        double sec;
        double hour;
        double day;
        double year;
        // distance
        double meter;
	double fm;
        double cm;
        double km;
        double fermi;
        double angstrom;
        double AU;
        double parsec;
        double kpc;
	// Area
	double barn;
        // luminocity
        double picobarn;
        double femtobarn;
        // presure
        double Pascal;
        double hPascal;
        double atm;
        double psi;
        // temperature
        double Kelvin;
        // angle
        double degree;
        double proton_mass;
        double neutron_mass;
	double Wboson_mass;
	double Zboson_mass;        
        PhysConst(void);
};

#endif
