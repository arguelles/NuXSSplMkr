#include "physconst.h"

PhysConst::PhysConst(void){
    /* PHYSICS CONSTANTS
    #===========================================================================
    # NAME
    #===========================================================================
    */

    name = "STD";                    // Default values
    linestyle = "solid";             // Default linestyle in plots
    markerstyle = "*";               // Default marker style
    colorstyle = "red";              // Default color style
    savefilename = "output.dat";     // Default color style

    /*
    #===============================================================================
    # MATH
    #===============================================================================
    */
    pi=3.14159265;		    // Pi
    piby2=1.5707963268;	            // Pi/2
    sqrt2=1.4142135624;	            // Sqrt[2]
    ln2 = log(2.0);                 // log[2]

    /*
    #===============================================================================
    # EARTH
    #===============================================================================
    */
    earthradius = 6371.0;	    // [km] Earth radius
    /*
    #===============================================================================
    # SUN
    #===============================================================================
    */
    sunradius = 109.0*earthradius;  // [km] Sun radius

    /*
    #===============================================================================
    # # PHYSICAL CONSTANTS
    #===============================================================================
    */
    GF = 1.16639e-23;	            // [eV^-2] Fermi Constant
    Na = 6.0221415e+23;		    // [mol cm^-3] Avogadro Number
    sw_sq = 0.2312;                 // [dimensionless] sin(th_weinberg) ^2
    G  = 6.67300e-11;               // [m^3 kg^-1 s^-2]
    alpha = 1.0/137.0;              // [dimensionless] fine-structure constant

    /*
    #===============================================================================
    # UNIT CONVERSION FACTORS
    #===============================================================================
    */
    // Energy
    TeV = 1.0e12;                   // [eV/TeV]
    GeV = 1.0e9;                    // [eV/GeV]
    MeV = 1.0e6;                    // [eV/MeV]
    keV = 1.0e3;                    // [eV/keV]
    Joule = 1/1.60225e-19;          // [eV/J]
    erg = (1.0e-7)*Joule;
    // Mass
    kg = 5.62e35;                   // [eV/kg]
    gr = 1e-3*kg;                   // [eV/g]
    // Time
    sec = 1.523e15;                 // [eV^-1/s]
    hour = 3600.0*sec;              // [eV^-1/h]
    day = 24.0*hour;                // [eV^-1/d]
    year = 365.0*day;               // [eV^-1/yr]
    // Distance
    meter = 5.076e6;                // [eV^-1/m]
    //fm = 1.0e-15*meter;		    // [ev^-1/fm]
    cm = 1.0e-2*meter;              // [eV^-1/cm]
    km = 1.0e3*meter;               // [eV^-1/km]
    fermi = 1.0e-15*meter;          // [eV^-1/fm]
    angstrom = 1.0e-10*meter;       // [eV^-1/A]
    AU = 149.60e9*meter;            // [eV^-1/AU]
    parsec = 3.08568025e16*meter;   // [eV^-1/parsec]
    kpc = 3.08568025e19*meter;   // [eV^-1/parsec]
    // Area
    barn = pow(10., -28)*meter*meter;   // [eV^-2/barn]
    // luminocity
    picobarn = 1.0e-36*pow(cm,2);       // [eV^-2/pb]
    femtobarn = 1.0e-39*pow(cm,2);      // [eV^-2/fb]
    // Presure
    Pascal = Joule/pow(meter,3);        // [eV^4/Pa]
    hPascal = 100.0*Pascal;         // [eV^4/hPa]
    atm = 101325.0*Pascal;          // [eV^4/atm]
    psi = 6893.0*Pascal;            // [eV^4/psi]
    // Temperature
    Kelvin = 1/1.1604505e4;         // [eV/K]
    // Angle
    degree = pi/180.0;              // [rad/degree]
    
    proton_mass = 0.938272*GeV;
    neutron_mass = 0.939565*GeV; 
    Wboson_mass = 80.385*GeV;
    Zboson_mass = 91.1876*GeV;

 
};

