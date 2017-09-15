#include "lhapdf_cross_section.h"

//#define DEBUG
//#define _DIPOLE__

LHAXS::LHAXS(std::string PDFname){
    //initialize Constants
    pdfname = PDFname;
    pc = new PhysConst();
    GF2 = SQ(pc->GF);
    M_iso  =    0.5*(pc->proton_mass + pc->neutron_mass);
    Mw2 = SQ(pc->Wboson_mass);
    Mz2 = SQ(pc->Zboson_mass);
    set = new LHAPDF::PDFSet(pdfname);
    if (set->errorType()=="custom")
        is_var = true;
    else
        is_var = false;
    
    nmem = set->size()-1;
    pdfs = set->mkPDFs();
    //error_band = 100*boost::math::erf(1./sqrt(2.));
    error_band = 95;
    //error_band = 68;
}

double LHAXS::SigR_Nu_LO_NC(double x, double y, map<int,LHAPDF::PDFUncertainty> dis, std::map<std::pair<int,int>,double> cov_m, int c){
  double k = 0;
	double q0 = 0.;
	double q0bar = 0.;

  //Following HEP PH 0407371 Eq. (21)
	double y_p = 1. + (1.- y) * (1. - y);
	double y_m = 1. - (1.- y) * (1. - y);
	double a = y_p + CP_factor*y_m;
	double b = y_p - CP_factor*y_m;

  double s_w = 0.2223;

  double Lu2 = ( 1. - (4./3.)*s_w) * ( 1. - (4./3.)*s_w);
  double Ld2 = (-1. + (2./3.)*s_w) * (-1. + (2./3.)*s_w);
  double Ru2 = (    - (4./3.)*s_w) * (    - (4./3.)*s_w);
  double Rd2 = (      (2./3.)*s_w) * (      (2./3.)*s_w);

  map<int,double> SigRcoef;
  if (CP_factor > 0 ){
    SigRcoef[1]  =    0.5*(Lu2 + Ld2) + 0.5*(Ru2 + Rd2)*(1. - y)*(1. - y);
    SigRcoef[-1] =    0.5*(Ru2 + Rd2) + 0.5*(Lu2 + Ld2)*(1. - y)*(1. - y);
    SigRcoef[2]  =    0.5*(Lu2 + Ld2) + 0.5*(Ru2 + Rd2)*(1. - y)*(1. - y);
    SigRcoef[-2] =    0.5*(Ru2 + Rd2) + 0.5*(Lu2 + Ld2)*(1. - y)*(1. - y);
    SigRcoef[3]  =    (Ld2 + Rd2) + (Ld2 + Rd2)*(1. - y)*(1. - y);
    SigRcoef[-3] =    (Ld2 + Rd2) + (Ld2 + Rd2)*(1. - y)*(1. - y);
    SigRcoef[4]  =    (Lu2 + Ru2) + (Lu2 + Ru2)*(1. - y)*(1. - y);
    SigRcoef[-4] =    (Lu2 + Ru2) + (Lu2 + Ru2)*(1. - y)*(1. - y);
    SigRcoef[5]  =    (Ld2 + Rd2) + (Ld2 + Rd2)*(1. - y)*(1. - y);
    SigRcoef[-5] =    (Ld2 + Rd2) + (Ld2 + Rd2)*(1. - y)*(1. - y);
    //SigRcoef[6]  =    (Lu2 + Ru2) + (Lu2 + Ru2)*(1. - y)*(1. - y);
    //SigRcoef[-6] =    (Lu2 + Ru2) + (Lu2 + Ru2)*(1. - y)*(1. - y);
    SigRcoef[21] =   0. ;
  } else {
    //fixes for antineutrinos
    SigRcoef[1]  =    0.5*(Lu2 + Ld2)*(1. - y)*(1. - y) + 0.5*(Ru2 + Rd2);
    SigRcoef[-1] =    0.5*(Ru2 + Rd2)*(1. - y)*(1. - y) + 0.5*(Lu2 + Ld2);
    SigRcoef[2]  =    0.5*(Lu2 + Ld2)*(1. - y)*(1. - y)+ 0.5*(Ru2 + Rd2);
    SigRcoef[-2] =    0.5*(Ru2 + Rd2)*(1. - y)*(1. - y)+ 0.5*(Lu2 + Ld2);
    SigRcoef[3]  =    (Ld2 + Rd2) *(1. - y)*(1. - y)+ (Ld2 + Rd2);
    SigRcoef[-3] =    (Ld2 + Rd2)*(1. - y)*(1. - y) + (Ld2 + Rd2);
    SigRcoef[4]  =    (Lu2 + Ru2)*(1. - y)*(1. - y) + (Lu2 + Ru2);
    SigRcoef[-4] =    (Lu2 + Ru2)*(1. - y)*(1. - y) + (Lu2 + Ru2);
    SigRcoef[5]  =    (Ld2 + Rd2)*(1. - y)*(1. - y) + (Ld2 + Rd2);
    SigRcoef[-5] =    (Ld2 + Rd2)*(1. - y)*(1. - y) + (Ld2 + Rd2);
    //SigRcoef[6]  =    (Lu2 + Ru2)*(1. - y)*(1. - y) + (Lu2 + Ru2);
    //SigRcoef[-6] =    (Lu2 + Ru2)*(1. - y)*(1. - y) + (Lu2 + Ru2);
    SigRcoef[21] =   0. ;
  }

	// mean value
    switch (c){
        case 0:
            for( int p : partons ) {
                k += SigRcoef[p]*dis[p].central;
                //std::cout << "p_pdf " << p << " " << dis[p].central << std::endl;
            }
            //std::cout << k << std::endl;
            break;

	// upper error
        case 1:
	        for( int p : partons ) {
	        	k += SQ(SigRcoef[p]*dis[p].errplus);
	        	for( int u : partons ) {
	        		if (dis[p].errplus*dis[u].errplus > 0.){
	        			k += ((SigRcoef[p]*SigRcoef[u])*
	        				(dis[p].errplus*dis[u].errplus)*
	        				cov_m[{p,u}]);
	        		}
	        	}
	        }
	        k = (1.0)*sqrt(k);
          for( int p : partons ) {
              k += SigRcoef[p]*dis[p].central;
          }
          break;

	// lower error
        case -1:
	        for( int p : partons ) {
	        	k += SQ(SigRcoef[p]*dis[p].errminus);
	        	for( int u : partons ) {
	        		if (dis[p].errminus*dis[u].errminus > 0.){
	        			k += ((SigRcoef[p]*SigRcoef[u])*
	        				(dis[p].errminus*dis[u].errminus)*
	        				cov_m[{p,u}]);
	        		}
	        	}
	        }
	        k = (-1.0)*sqrt(k);
          for( int p : partons ) {
              k += SigRcoef[p]*dis[p].central;
          }
          break;

    // symmetric
        case 9:
	        for( int p : partons ) {
	        	k += SQ(SigRcoef[p]*dis[p].errsymm);
	        	for( int u : partons ) {
	        		if (dis[p].errsymm*dis[u].errsymm > 0.){
	        			k += ((SigRcoef[p]*SigRcoef[u])*
	        				(dis[p].errsymm*dis[u].errsymm)*
	        				cov_m[{p,u}]);
	        		}
	        	}
	        }	
	        k = sqrt(k);
            break;
    }

	return k;
}

double LHAXS::SigR_Nu_LO(double x, double y, map<int,LHAPDF::PDFUncertainty> dis, std::map<std::pair<int,int>,double> cov_m, int c){
	double k = 0.;

  //Following HEP PH 0407371 Eq. (21)
	double y_p = 1. + (1.- y) * (1. - y);
	double y_m = 1. - (1.- y) * (1. - y);
	//double y_p = (1. - d_lepton / x) + (1.- d_lepton/x - y) * (1. - y);
	//double y_m = (1. - d_lepton / x) - (1.- d_lepton/x - y) * (1. - y);
	double a = y_p + CP_factor*y_m;
	double b = y_p - CP_factor*y_m;
  //std::cout << "x " << x << " y " << y << std::endl;
  //std::cout << "d_lepton " << d_lepton << " CP_factor " << CP_factor<< std::endl;
  //std::cout << "a " << a << " b " << b << std::endl;

  map<int,double> SigRcoef;
	SigRcoef[1]  =    a ;
	SigRcoef[-1] =    b ;
	SigRcoef[2]  =    a ;
	SigRcoef[-2] =    b ;
	SigRcoef[3]  = 2.*a ;
	SigRcoef[-3] =   0. ;
	SigRcoef[4]  =   0. ;
	SigRcoef[-4] = 2.*b ;
	SigRcoef[5]  = 2.*a ;
	SigRcoef[-5] =   0. ;
	SigRcoef[21] =   0. ;

  if (CP_factor < 0 ){
    //fixes for antineutrinos
    SigRcoef[3]   = 0. ;
    SigRcoef[-3]  = 2.*b ;
    SigRcoef[4]   = 2.*a ;
    SigRcoef[-4]  = 0. ;
    SigRcoef[5]   = 0. ;
    SigRcoef[-5]  = 2.*b ;
  }
//  std::cout << "x " << x << " y " << y <<" a " << a << " b " << b << " CPF " << CP_factor <<std::endl;

	// mean value
    switch (c){
        case 0:
            for( int p : partons ) {
                k += SigRcoef[p]*dis[p].central;
                //std::cout << "Sigrcoef " << SigRcoef[p] << " p_pdf " << p << " " << dis[p].central << std::endl;
            }
            //std::cout << k << std::endl;
            break;

	// upper error
        case 1:
	        for( int p : partons ) {
	        	k += SQ(SigRcoef[p]*dis[p].errplus);
	        	for( int u : partons ) {
	        		if (dis[p].errplus*dis[u].errplus > 0.){
	        			k += ((SigRcoef[p]*SigRcoef[u])*
	        				(dis[p].errplus*dis[u].errplus)*
	        				cov_m[{p,u}]);
	        		}
	        	}
	        }
	        k = (1.0)*sqrt(k);
          for( int p : partons ) {
              k += SigRcoef[p]*dis[p].central;
          }
          break;

	// lower error
        case -1:
	        for( int p : partons ) {
	        	k += SQ(SigRcoef[p]*dis[p].errminus);
	        	for( int u : partons ) {
	        		if (dis[p].errminus*dis[u].errminus > 0.){
	        			k += ((SigRcoef[p]*SigRcoef[u])*
	        				(dis[p].errminus*dis[u].errminus)*
	        				cov_m[{p,u}]);
	        		}
	        	}
	        }
	        k = (-1.0)*sqrt(k);
          for( int p : partons ) {
              k += SigRcoef[p]*dis[p].central;
          }
          break;

    // symmetric
        case 9:
	        for( int p : partons ) {
	        	k += SQ(SigRcoef[p]*dis[p].errsymm);
	        	for( int u : partons ) {
	        		if (dis[p].errsymm*dis[u].errsymm > 0.){
	        			k += ((SigRcoef[p]*SigRcoef[u])*
	        				(dis[p].errsymm*dis[u].errsymm)*
	        				cov_m[{p,u}]);
	        		}
	        	}
	        }	
	        k = sqrt(k);
            break;
    }

	return k;
}

//==================================================================================
// structure functions
//==================================================================================

double LHAXS::F1(double x, double q2){
    // only true at LO
    return F2(x,q2)/(2.*x);
}

double LHAXS::F2(double x, double q2){
    // only true at LO
    auto s = PDFExtract(x,q2);
    return F2(s);
}

double LHAXS::F3(double x, double q2){
    return xF3(x,q2)/x;
}
double LHAXS::xF3(double x, double q2){
    // only true at LO
    auto s = PDFExtract(x,q2);
    return xF3(s);
}

double LHAXS::F4(double x, double q2){
    // valid at all orders
    //return (F2(x,q2)/(2.x) - F1(x,q2))/2.0;
    // only true at LO
    return 0.;
}

double LHAXS::F5(double x, double q2){
    // valid at all orders
    return F2(x,q2)/(2.*x);
}

double LHAXS::F2(map<int, double>& xq_arr) const{
    double k = 0.;

	map<int,double> F2coef;
	F2coef[1]  = 1.;
	F2coef[-1] = 1.;
	F2coef[2]  = 1.;
	F2coef[-2] = 1.;
	F2coef[3]  = 2.;
	F2coef[-3] = 0.;
	F2coef[4]  = 0.;
	F2coef[-4] = 2.;
	F2coef[5]  = 2.;
	F2coef[-5] = 0.;
	F2coef[21] = 0.;

	// mean value
	for( int p : partons ) {
		k += F2coef[p]*xq_arr[p];
	}	
    return k;
}

/*
double LHAXS::F2NLO(const map<int,double>* xq_arr) const {

}
*/

double LHAXS::xF3(map<int, double>& xq_arr) const{
	double k=0.;

	// xF3 coeficients
	map<int,double> F3coef;
	F3coef[1]  = 1.;
	F3coef[-1] = -1.;
	F3coef[2]  = 1.;
	F3coef[-2] = -1.;
	F3coef[3]  = 2.;
	F3coef[-3] = 0.;
	F3coef[4]  = 0.;
	F3coef[-4] = -2.;
	F3coef[5]  = 2.;
	F3coef[-5] = 0.;
	F3coef[21] = 0.;

	// mean value
	for( int p : partons ) {
		k += F3coef[p]*xq_arr[p];
	}	

    return k;
}

//==================================================================================
// reduced cross section
//==================================================================================

map<int,double> LHAXS::PDFExtract(double x, double q2){
    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(pdfs[0]);
    string xt = "nearest";
    grid_central -> setExtrapolator(xt);

    map<int,double> xq_arr;
    for ( int p : partons ){
      xq_arr[p] = grid_central -> xfxQ2(p,x,q2/SQ(pc->GeV));
    }
    return xq_arr;
}

double LHAXS::SigRed_Evaluate(double q2, double x, double y){
    //Take this representation of q !!
    double y_p = (1. - d_lepton / x) + (1.- d_lepton/x - y) * (1. - y);
    double y_m = (1. - d_lepton / x) - (1.- d_lepton/x - y) * (1. - y);

    map<int,double> xqa = PDFExtract(x,q2);

    return y_p*F2(xqa) + CP_factor*y_m*xF3(xqa);
}

double LHAXS::SigRed_TMC(double x, double y, double q2){
    double norm_tmc = GF2*M_iso*ENU/(M_PI*SQ(1.+q2/M_boson2));

    double sum_terms = 0.0;

    sum_terms += (y*y*x + SQ(M_lepton)*y / (2.*ENU*M_iso)) * F1_TMC(x,q2);
    sum_terms += (1. - SQ(M_lepton) / (4.*SQ(ENU)) - y + M_iso*x*y / (2.*ENU)) * F2_TMC(x,q2);
    sum_terms += (x*y - x*y*y / 2. - SQ(M_lepton)*y / (4.*ENU*M_iso)) * F3_TMC(x,q2) * CP_factor;
    sum_terms += (SQ(M_lepton)*(SQ(M_lepton) + q2) / (4.*SQ(ENU*M_iso)*x) ) * F4_TMC(x,q2);
    sum_terms += -(SQ(M_lepton) / (ENU*M_iso) ) * F5_TMC(x,q2);

    //std::cout << norm_tmc << " "  << ENU << " "  << M_iso << " "  << M_boson2 << " "  << F1_TMC(x, q2) <<  " " <<Mw2 << std::endl;

    return norm_tmc * sum_terms;
}

double LHAXS::SigR_Nu_LO_NC(double x,double y, map<int, double> xq_arr){
  double s_w = 0.2223;

  double Lu2 = ( 1. - (4./3.)*s_w) * ( 1. - (4./3.)*s_w);
  double Ld2 = (-1. + (2./3.)*s_w) * (-1. + (2./3.)*s_w);
  double Ru2 = (    - (4./3.)*s_w) * (    - (4./3.)*s_w);
  double Rd2 = (      (2./3.)*s_w) * (      (2./3.)*s_w);

  double u = xq_arr[1];
  double d = xq_arr[2];
  double s = xq_arr[3];
  double c = xq_arr[4];
  double b = xq_arr[5];
  //double t = xq_arr[6];
  double ubar = xq_arr[-1];
  double dbar = xq_arr[-2];
  double sbar = xq_arr[-3];
  double cbar = xq_arr[-4];
  double bbar = xq_arr[-5];
  //double tbar = xq_arr[-6];

  //double q0   = 0.5*(u + d)*(Lu2 + Ld2) + 0.5*(ubar + dbar)*(Ru2 + Rd2) + (s + b + sbar + bbar)*(Ld2 + Rd2) + (cbar + c + tbar + t)*(Lu2 + Ru2);
  //double q0bar= 0.5*(u + d)*(Ru2 + Rd2) + 0.5*(ubar + dbar)*(Lu2 + Ld2) + (s + b + sbar + bbar)*(Ld2 + Rd2) + (cbar + c + tbar + t)*(Lu2 + Ru2);

  double q0   = 0.5*(u + d)*(Lu2 + Ld2) + 0.5*(ubar + dbar)*(Ru2 + Rd2) + (s + b + sbar + bbar)*(Ld2 + Rd2) + (cbar + c)*(Lu2 + Ru2);
  double q0bar= 0.5*(u + d)*(Ru2 + Rd2) + 0.5*(ubar + dbar)*(Lu2 + Ld2) + (s + b + sbar + bbar)*(Ld2 + Rd2) + (cbar + c)*(Lu2 + Ru2);

  if (CP_factor > 0 )
    return (q0 + q0bar*(1.-y)*(1.-y));
  else
    return (q0*(1.-y)*(1.-y)+ q0bar);
}

double LHAXS::SigR_Nu_LO(double x, double y, map<int,double> xq_arr){
	double k = 0.;

	double y_p = (1. - d_lepton / x) + (1.- d_lepton/x - y) * (1. - y);
	double y_m = (1. - d_lepton / x) - (1.- d_lepton/x - y) * (1. - y);
	double a = y_p + CP_factor*y_m;
	double b = y_p - CP_factor*y_m;

  map<int,double> SigRcoef;

  // Coefficients for CC
	SigRcoef[1]  =    a ;
	SigRcoef[-1] =    b ;
	SigRcoef[2]  =    a ;
	SigRcoef[-2] =    b ;
	SigRcoef[3]  = 2.*a ;
	SigRcoef[-3] =   0. ;
	SigRcoef[4]  =   0. ;
	SigRcoef[-4] = 2.*b ;
	SigRcoef[5]  = 2.*a ;
	SigRcoef[-5] =   0. ;
	SigRcoef[21] =   0. ;

  if (CP_factor < 0 ){
    //fixes for antineutrinos
    SigRcoef[3]   = 0. ;
    SigRcoef[-3]  = 2.*b ;
    SigRcoef[4]   = 2.*a ;
    SigRcoef[-4]  = 0. ;
    SigRcoef[5]   = 0. ;
    SigRcoef[-5]  = 2.*b ;
  }

  for( int p : partons ) {
      k += SigRcoef[p]*xq_arr[p];
      //std::cout << "p_pdf " << p << " " << xq_arr[p] << std::endl;
  }
  //std::cout << k << std::endl;

  return k;
}

double LHAXS::EvaluateVar(double Q2, double x, double y, int var){
    // only evaluates central values
    double q = sqrt(Q2)/pc->GeV;

    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(pdfs[var]);
    string xt = "nearest";
    grid_central -> setExtrapolator(xt);

    map<int,double> xq_arr;
    for (int p : partons){
      xq_arr[p] = grid_central -> xfxQ(p,x,q);
    }

    if(INT_TYPE==CC)
      return SigR_Nu_LO(x, y, xq_arr);
    else
      return SigR_Nu_LO_NC(x, y, xq_arr);
}

double LHAXS::Evaluate(double Q2, double x, double y){
    // only evaluates central values
    double q = sqrt(Q2)/pc->GeV;

    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(pdfs[0]);
    string xt = "nearest";
    grid_central -> setExtrapolator(xt);

    map<int,double> xq_arr;
    for ( int p : partons ){
      xq_arr[p] = grid_central -> xfxQ(p,x,q);
    }

    if(INT_TYPE==CC)
      return SigR_Nu_LO(x, y, xq_arr);
    else
      return SigR_Nu_LO_NC(x, y, xq_arr);
}

double LHAXS::Evaluate (double Q2, double x, double y, int a){
    double q = sqrt(Q2)/pc->GeV;
    map<int,vector<double>> xpdf_arr;
    map<int,LHAPDF::PDFUncertainty> xer_arr;
    int gluon_index = 21;
    //std::cout << "q_gev " << q << std::endl;

    //Quark Flavors
    for(int p : partons){
        if (p==gluon_index)
            continue;

        vector<double> xf;
        vector<double> delta_xf;

        LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(pdfs[0]);
        string xt = "nearest";
        grid_central -> setExtrapolator(xt);

        for (size_t imem=0; imem<=nmem; imem++){
            LHAPDF::GridPDF* grid = dynamic_cast<LHAPDF::GridPDF*>(pdfs[imem]);
            grid -> setExtrapolator(xt);
            xf.push_back(grid -> xfxQ(p, x, q));
            double hij = grid -> xfxQ(p, x, q) - grid_central -> xfxQ(p,x,q);
            delta_xf.push_back(hij);
        }

        LHAPDF::PDFUncertainty xerror;
        if ( is_var ){
            double smp = 0.;
            double smm = 0.;
            double spp,spm;

            // model error
            for ( int ij = 1 ; ij < 9; ij ++) {
                double dfx = delta_xf[ij];
                if (dfx > 0.)
                    smp += dfx*dfx;
                else
                    smm += dfx*dfx;
            }	

            // parametrization error
            double dfxpm = 0;
            double dfxmm = 0;
            for ( int ij = 9 ; ij < delta_xf.size(); ij ++) {
                double dfx = delta_xf[ij];
                if (dfx > 0. && dfx > dfxpm)
                {
                    spp = dfx*dfx;
                    dfxpm = dfx;
                }
                else if ( dfx < 0. && dfx < dfxmm)
                {
                    spm = dfx*dfx;
                    dfxmm = dfx;
                }
            }	

            // add
            xerror.central = grid_central -> xfxQ(p,x,q);
            xerror.errplus = sqrt(smp+spp);
            xerror.errminus = sqrt(smm+spm);
            xerror.errsymm = max(xerror.errplus,xerror.errminus);
        } else 
            xerror = set->uncertainty(xf,error_band);

		xer_arr.insert({p,xerror});
		xpdf_arr.insert({p,xf});
		parton_num[p] = xer_arr.size()-1;
    } // end of quark loop

    // Gluon
    vector<double> xg;
    vector<double> delta_xg;
    string xt = "nearest";
    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(pdfs[0]);
    grid_central -> setExtrapolator(xt);

    for (size_t imem=0; imem<=nmem; imem++){
        LHAPDF::GridPDF* grid = dynamic_cast<LHAPDF::GridPDF*>(pdfs[imem]);
        grid -> setExtrapolator(xt);
        double hij =  pdfs[imem] -> xfxQ(gluon_index, x, q) - pdfs[0] -> xfxQ(gluon_index,x,q);
        delta_xg.push_back(hij);
        xg.push_back(pdfs[imem] -> xfxQ(gluon_index, x, q));
    }

    LHAPDF::PDFUncertainty xerror;

    if ( is_var ){
        double smp = 0.;
        double smm = 0.;
        double spp,spm;
        // model error
        for ( int ij = 1 ; ij < 9; ij ++) {
            double dfx = delta_xg[ij];
            if (dfx > 0.)
                smp += dfx*dfx;
            else
                smm += dfx*dfx;
        }	
        // parametrization error
        double dfxpm = 0;
        double dfxmm = 0;
        for ( int ij = 9 ; ij < delta_xg.size(); ij ++) {
            double dfx = delta_xg[ij];
            if (dfx > 0. && dfx > dfxpm)
            {
                spp = dfx*dfx;
                dfxpm = dfx;
            }
            else if ( dfx < 0. && dfx < dfxmm)
            {
                spm = dfx*dfx;
                dfxmm = dfx;
            }
        }	
        xerror.central = grid_central -> xfxQ(21,x,q);
        xerror.errplus = sqrt(smp+spp);
        xerror.errminus = sqrt(smm+spm);
        xerror.errsymm = max(xerror.errplus,xerror.errminus);
    } else 
        xerror = set->uncertainty(xg,error_band);

	xpdf_arr.insert({gluon_index,xg});
	xer_arr.insert({gluon_index,xerror});
	parton_num[gluon_index] = xer_arr.size()-1;

  std::map<std::pair<int,int>,double> cor_mat;
  for(int kk : partons){
    for(int jj : partons){
      if (jj == kk || is_var)
        cor_mat.insert({{kk,jj},0.0});
      else
        cor_mat.insert({{kk,jj},set->correlation(xpdf_arr[kk],xpdf_arr[jj])});
    }
  }

  //std::cout << xpdf_arr.size() << std::endl;
  //std::cout << xpdf_arr[0].size() << std::endl;
  /*
  for( int p : partons){
     std::cout << " ppdf " << p << " " << xpdf_arr[p][0] << " " << xer_arr[p].central;
  }
  std::cout << std::endl;
  */
  if(INT_TYPE==CC)
    return SigR_Nu_LO(x, y, xer_arr, cor_mat, a);
  else
    return SigR_Nu_LO_NC(x, y, xer_arr, cor_mat, a);
}

/*
double LHAXS::Evaluate (double Q2, double x, double y, int a){
    double q = sqrt(Q2/pc->GeV);
    vector<vector<double>> xpdf_arr;
    vector<LHAPDF::PDFUncertainty> xer_arr;

    //Quark Flavors
    for (int i = -5; i<=5;i++){
        vector<double> xf;
        vector<double> delta_xf;
        if (i==0)
            continue;

        LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(pdfs[0]);
        string xt = "nearest";
        grid_central -> setExtrapolator(xt);

        for (size_t imem=0; imem<=nmem; imem++){

            LHAPDF::GridPDF* grid = dynamic_cast<LHAPDF::GridPDF*>(pdfs[imem]);
            grid -> setExtrapolator(xt);
            xf.push_back(grid -> xfxQ(i, x, q));
            double hij = + grid -> xfxQ(i, x, q) - grid_central -> xfxQ(i,x,q);
            delta_xf.push_back(hij);
        }

        LHAPDF::PDFUncertainty xerror;
        if ( is_var ){
            double smp = 0.;
            double smm = 0.;
            double spp,spm;

            // model error
            for ( int ij = 1 ; ij < 9; ij ++) {
                double dfx = delta_xf[ij];
                if (dfx > 0.)
                    smp += dfx*dfx;
                else
                    smm += dfx*dfx;
            }	

            // parametrization error
            double dfxpm = 0;
            double dfxmm = 0;
            for ( int ij = 9 ; ij < delta_xf.size(); ij ++) {
                double dfx = delta_xf[ij];
                if (dfx > 0. && dfx > dfxpm)
                {
                    spp = dfx*dfx;
                    dfxpm = dfx;
                }
                else if ( dfx < 0. && dfx < dfxmm)
                {
                    spm = dfx*dfx;
                    dfxmm = dfx;
                }
            }	

            // add
            xerror.central = grid_central -> xfxQ(i,x,q);
            xerror.errplus = sqrt(smp+spp);
            xerror.errminus = sqrt(smm+spm);
            xerror.errsymm = max(xerror.errplus,xerror.errminus);
        } else 
            xerror = set->uncertainty(xf,error_band);

		xer_arr.push_back(xerror);
		xpdf_arr.push_back(xf);
		parton_num[i] = xer_arr.size()-1;
    } // end of quark loop

    // Gluon
    vector<double> xg;
    vector<double> delta_xg;
    string xt = "nearest";
    LHAPDF::GridPDF* grid_central = dynamic_cast<LHAPDF::GridPDF*>(pdfs[0]);
    grid_central -> setExtrapolator(xt);

    for (size_t imem=0; imem<=nmem; imem++){
        LHAPDF::GridPDF* grid = dynamic_cast<LHAPDF::GridPDF*>(pdfs[imem]);
        grid -> setExtrapolator(xt);
        double hij =  pdfs[imem] -> xfxQ(21, x, q) - pdfs[0] -> xfxQ(21,x,q);
        delta_xg.push_back(hij);
        xg.push_back(pdfs[imem] -> xfxQ(21, x, q));
    }

    LHAPDF::PDFUncertainty xerror;

    if ( is_var ){
        double smp = 0.;
        double smm = 0.;
        double spp,spm;
        // model error
        for ( int ij = 1 ; ij < 9; ij ++) {
            double dfx = delta_xg[ij];
            if (dfx > 0.)
                smp += dfx*dfx;
            else
                smm += dfx*dfx;
        }	
        // parametrization error
        double dfxpm = 0;
        double dfxmm = 0;
        for ( int ij = 9 ; ij < delta_xg.size(); ij ++) {
            double dfx = delta_xg[ij];
            if (dfx > 0. && dfx > dfxpm)
            {
                spp = dfx*dfx;
                dfxpm = dfx;
            }
            else if ( dfx < 0. && dfx < dfxmm)
            {
                spm = dfx*dfx;
                dfxmm = dfx;
            }
        }	
        xerror.central = grid_central -> xfxQ(21,x,q);
        xerror.errplus = sqrt(smp+spp);
        xerror.errminus = sqrt(smm+spm);
        xerror.errsymm = max(xerror.errplus,xerror.errminus);
    } else 
        xerror = set->uncertainty(xg,error_band);

	xer_arr.push_back(xerror);
	xpdf_arr.push_back(xg);
	parton_num[21] = xer_arr.size()-1;

    vector<vector<double>> cor_mt;
    //if (a!==0){
	    for(int kk = 0 ; kk < xpdf_arr.size(); kk++){
	    	vector<double> cor_v;
	    	for(int jj = 0 ; jj < xpdf_arr.size(); jj++){
	    		if (jj == kk || is_var)
	    			cor_v.push_back(0.0);
	    		else
	    			cor_v.push_back(set->correlation(xpdf_arr[kk],xpdf_arr[jj]));
	    	}
	    	cor_mt.push_back(cor_v);
	    }
   //}
      std::cout << xpdf_arr.size() << std::endl;
      std::cout << xpdf_arr[0].size() << std::endl;
   for( int p : partons){
     std::cout << " ppdf " << p << " " << xpdf_arr[p][0];
   }
   std::cout << std::endl;
    return SigR_Nu_LO(x, y, xer_arr, cor_mt, a);
}
*/

//==================================================================================
// SET OPTIONS 
//==================================================================================

void LHAXS::Set_QCDOrder(QCDOrder qcdorder_){
  qcdorder = qcdorder_;
}

void LHAXS::Set_M_Lepton(double ml){
    M_lepton = ml;
}

void LHAXS::Set_CP_factor(double cp){
    CP_factor = cp;
}

void LHAXS::Set_InteractionType(Current c){
    INT_TYPE = c;
    if( INT_TYPE == CC)
        M_boson2 = Mw2;
    else
        M_boson2 = Mz2;
}

void LHAXS::Set_Neutrino_Energy(double enu){
  ENU = enu;
  ienu = true;
}

void LHAXS::Set_Variant(int var){
  ivar=var;
}

//==================================================================================
// TMC CORRECTIONS
// HEP-PH 0307023v2 (S. Kretzer, M.H. Reno)
// HEP-PH 07091775
//==================================================================================

double LHAXS::R(double x, double q2){
    return sqrt(1. + 4.*x*x*SQ(M_iso)/q2);
}

double LHAXS::Xi(double x, double q2){
    double denumerator = 1. + R(x,q2);
    return 2.*x / denumerator;
}


template<class T,double (T::*f)(double,double),int n,int m>
double LHAXS::HGeneric(double xi, double q2){
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    double result, error;
    
    gsl_function F;
    F.function = &HK<T,f,n,m>;
    F.params = this;
    // set q2
    Q2 = q2;
    gsl_integration_qags ( &F, xi, 1, 0, 1.e-7, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

double LHAXS::H1(double xi, double q2){
   return HGeneric<LHAXS,&LHAXS::F1,1,2>(xi,q2);
}

double LHAXS::H2(double xi, double q2){
   return HGeneric<LHAXS,&LHAXS::F2,2,1>(xi,q2);
}

double LHAXS::H3(double xi, double q2){
   return HGeneric<LHAXS,&LHAXS::F3,1,1>(xi,q2);
}

double LHAXS::H4(double xi, double q2){
   return HGeneric<LHAXS,&LHAXS::F4,1,4>(xi,q2);
}

double LHAXS::H5(double xi, double q2){
   return HGeneric<LHAXS,&LHAXS::F5,1,2>(xi,q2);
}

double LHAXS::G2(double xi, double q2){
   return HGeneric<LHAXS,&LHAXS::F2,1,1>(xi,q2) - xi*HGeneric<LHAXS,&LHAXS::F2,2,1>(xi,q2);
}

double LHAXS::F1_TMC(double x, double q2){
    double xi = Xi(x, q2);
    double r  = R(x, q2);
    double g2 = G2(xi, q2);
    double h2 = H2(xi, q2);

    double sum1 = x/(xi * r) * F1(xi, q2);
    double sum2 = M_iso*M_iso* x*x  * h2 / (q2* r*r);
    double sum3 = 2.*SQ(M_iso*M_iso)* x*SQ(x) * g2 / (SQ(q2) * r*r*r );

    return sum1 + sum2 + sum3;
}

double LHAXS::F2_TMC(double x, double q2){
    double xi = Xi(x, q2);
    double r  = R(x, q2);
    double g2 = G2(xi, q2);
    double h2 = H2(xi, q2);
    
    double sum1 = x*x/(SQ(xi) * r*r*r) * F2(xi, q2);
    double sum2 = 6.*M_iso*M_iso* x*x*x  * h2 / (q2* r*r*r*r);
    
    double sum3 = 12.*M_iso*M_iso* SQ(SQ(x)) * g2 / (SQ(q2) * r*r*r*r*r );
    return sum1+sum2+sum3;
}

double LHAXS::F3_TMC(double x, double q2){
    double xi = Xi(x, q2);
    double r  = R(x, q2);
    double g2 = G2(xi, q2);
    double h3 = H3(xi, q2);
    
    double sum1 = x/(xi * r*r) * F3(xi, q2);
    double sum2 = 2.*M_iso*M_iso* x*x  * h3 / (q2* r*r*r);

    return sum1+sum2;
}

double LHAXS::xF3_TMC(double x, double q2){
    return x*F3_TMC(x,q2);
}
double LHAXS::F4_TMC(double x, double q2){
    return (F2_TMC(x,q2) / (4.*x)) - (F1_TMC(x,q2) / 2.);
}

double LHAXS::F5_TMC(double x, double q2){
    return F2_TMC(x,q2) / (2.*x);
}

double LHAXS::FL_TMC(double x, double q2){
    double r  = R(x, q2);
    return r*r*F2_TMC(x,q2) - 2.*x*F1_TMC(x,q2);
}
// helper functions

double LHAXS::KernelXS(double * k,int a){
  if (!ienu)
    throw std::runtime_error("energy not initialize");
  double x = k[0];
  double y = k[1];
  double Q2 = (2.*M_iso*ENU + SQ(M_iso))*x*y;

  // same for CC and NC
  double denum    = SQ(1. + Q2/M_boson2);
  double norm     = GF2*M_iso*ENU/(2.*M_PI*denum);

  //std::cout << ENU << " " << Q2  << " " <<  norm << " " << M_iso << " " << denum << " " << M_boson2 << std::endl;
  //Following HEP PH 0407371 Eq. (7)
  //if(INT_TYPE==NC)
  //  d_lepton=0;

  if(INT_TYPE==CC){
    double h = x*y + d_lepton;
    if ((1. + x*d_nucleon) * h*h - (x+ d_lepton)*h + x * d_lepton > 0.)
        return 0.;
  }
  return norm*Evaluate(Q2, x, y, a);
}

double LHAXS::KernelXSVar(double * k){
  double x = k[0];
  double y = k[1];
  double Q2 = (2.*M_iso*ENU + SQ(M_iso))*x*y;

  double denum    = SQ(1. + Q2/M_boson2);
  double norm     = GF2*M_iso*ENU/(2.*M_PI*denum);
  //std::cout << Evaluate(Q2, x, y, 0) << std::endl;

  d_lepton = SQ(M_lepton)/(2.*M_iso*ENU);

  if(INT_TYPE==CC){
    //Following HEP PH 0407371 Eq. (7)
    double h = x*y + d_lepton;
    if ((1. + x* d_nucleon) * h*h - (x+ d_lepton)*h + x * d_lepton > 0.)
        return 0.;
  }
  if(ivar==0)//central set
    return norm*Evaluate(Q2, x, y);
  else
    return norm*Evaluate(Q2, x, y, ivar);
}

double LHAXS::KernelXS(double * k){
  double x = k[0];
  double y = k[1];
  double Q2 = (2.*M_iso*ENU + SQ(M_iso))*x*y;

  double denum    = SQ(1. + Q2/M_boson2);
  double norm     = GF2*M_iso*ENU/(2.*M_PI*denum);
  //std::cout << Evaluate(Q2, x, y, 0) << std::endl;

  d_lepton = SQ(M_lepton)/(2.*M_iso*ENU);

  //Following HEP PH 0407371 Eq. (7)
  if(INT_TYPE==CC){
    double h = x*y + d_lepton;
    if ((1. + x* d_nucleon) * h*h - (x+ d_lepton)*h + x * d_lepton > 0.)
        return 0.;
  }
  return norm*Evaluate(Q2, x, y);
}

double LHAXS::KernelXS_TMC(double * k){
  double x = k[0];
  double y = k[1];
  double q2 = (2.*M_iso*ENU + SQ(M_iso))*x*y;

  //Following HEP PH 0407371 Eq. (7)
  double h = x*y + d_lepton;
  double cc = (1. + x* d_nucleon) * h*h - (x+ d_lepton)*h + x * d_lepton;
  //std::cout << x << " " << y << " " << q2 << "\t" << h << " " << cc << " ==> " << SigRed_TMC(x, y, q2) << std::endl;
  if ( cc > 0.)
      return 0.;
  return SigRed_TMC(x,y,q2);
}

double LHAXS::KernelXS_dsdyVar(double x){
    double s = 2.*M_iso*ENU + SQ(M_iso);
    //cout << s << " " << x << " " << Y_EMU << endl;
    double q2 = s*x*Y;

    double denum = SQ(1. + q2/M_boson2);
    double norm = GF2*M_iso*ENU/(2.*M_PI*denum);
    d_lepton = SQ(M_lepton)/(2.*M_iso*ENU);

    //std::cout << norm * Evaluate(q2, x, Y, ivar) << std::endl;
    //std::cout << Evaluate(q2, x, Y, ivar) << std::endl;
    return norm * Evaluate(q2, x, Y, ivar);
}

double LHAXS::KernelXS_dsdy(double x){
    double s = 2.*M_iso*ENU + SQ(M_iso);
    //cout << s << " " << x << " " << Y_EMU << endl;
    double q2 = s*x*Y;

    double denum = SQ(1. + q2/M_boson2);
    double norm = GF2*M_iso*ENU/(2.*M_PI*denum);

    d_lepton = SQ(M_lepton)/(2.*M_iso*ENU);

    return norm * Evaluate (q2, x, Y);
}

double LHAXS::dsdyVar(double y){
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
    double result, error;

    gsl_function F;
    F.function = &KernelHelper<LHAXS,&LHAXS::KernelXS_dsdyVar>;
    F.params = this;
    // set y
    Y = y;
    gsl_integration_qag ( &F, 0, 1, 0, 1.e-5, 5000, 6, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

double LHAXS::dsdy(double y){
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
    double result, error;

    gsl_function F;
    F.function = &KernelHelper<LHAXS,&LHAXS::KernelXS_dsdy>;
    F.params = this;
    // set y
    Y = y;
    gsl_integration_qag ( &F, 0, 1, 0, 1.e-5, 5000, 6, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

double LHAXS::totalVar(){
  double res,err;
  const unsigned long dim = 2; int calls = 50000;
  double xl[dim] = {  0.0 , 0.0 };
  double xu[dim] = {  1.0 , 1.0 };

  gsl_rng_env_setup ();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc (T);

  gsl_monte_function F = { &KernelHelper<LHAXS,&LHAXS::KernelXSVar>, dim, this};
  gsl_monte_vegas_state *s_vegas = gsl_monte_vegas_alloc (dim);

  // training
  //std::cout << "s_vegas: " << s_vegas << std::endl;
  //std::cout << &xl << " " << &xu << " "  << dim << " "  << calls << " "  << &r << std::endl;
  gsl_monte_vegas_integrate (&F, xl, xu, dim, 10000, r, s_vegas,
                              &res, &err);
  //std::cout << "stop" << std::endl;
  do
  {
  std::cout << "Here: "<<gsl_monte_vegas_chisq (s_vegas) -1. << std::endl;
  gsl_monte_vegas_integrate (&F, xl, xu, dim, calls, r, s_vegas,
                              &res, &err);
  }
  while (fabs (gsl_monte_vegas_chisq (s_vegas) - 1.0) > 0.5 );

  gsl_monte_vegas_free (s_vegas);
  gsl_rng_free (r);
  //std::cout << "Result: " << res << std::endl;
  return res;
}
double LHAXS::total(){
  double res,err;
  const unsigned long dim = 2; int calls = 50000;
  double xl[dim] = {  0.0 , 0.0 };
  double xu[dim] = {  1.0 , 1.0 };

  gsl_rng_env_setup ();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc (T);

  gsl_monte_function F = { &KernelHelper<LHAXS,&LHAXS::KernelXS>, dim, this};
  gsl_monte_vegas_state *s_vegas = gsl_monte_vegas_alloc (dim);

  // training
  //std::cout << "s_vegas: " << s_vegas << std::endl;
  //std::cout << &xl << " " << &xu << " "  << dim << " "  << calls << " "  << &r << std::endl;
  gsl_monte_vegas_integrate (&F, xl, xu, dim, 10000, r, s_vegas,
                              &res, &err);
  //std::cout << "stop" << std::endl;
  do
  {
  //std::cout << "Here: "<<gsl_monte_vegas_chisq (s_vegas) << std::endl;
  gsl_monte_vegas_integrate (&F, xl, xu, dim, calls, r, s_vegas,
                              &res, &err);
  }
  while (fabs (gsl_monte_vegas_chisq (s_vegas) - 1.0) > 0.5 );

  gsl_monte_vegas_free (s_vegas);
  gsl_rng_free (r);
  //std::cout << "Result: " << res << std::endl;
  return res;
}
