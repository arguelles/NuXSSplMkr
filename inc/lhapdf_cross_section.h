#ifndef __LHAPDF_XS_
#define __LHAPDF_XS_

#define SQ(X)  ((X)*(X))
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
#include "LHAPDF/Extrapolator.h"
#include "physconst.h"
//#include <boost/random.hpp>
#include "tools.h"
#include <vector>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <functional>
#include <cmath>
#include <map>

enum QCDOrder {LO,NLO,NNLO};
enum Current {CC,NC};

template<class T,double (T::*f)(double*)>
double KernelHelper(double* x,size_t dim, void* param){
  T* p = (T*) param;
  return (p->*f)(x);
}

template<class T,double (T::*f)(double)>
double KernelHelper(double x, void* param){
  T* p = (T*) param;
  return (p->*f)(x);
}

template<class T,double (T::*f)(double,double),int n,int m>
double HK(double x, void* param){
  T* p = (T*) param;
  return ((double)m)*((p->*f)(x,p->Q2))/pow(x,(double)n);
}

class LHAXS{
    private:
        double s_w,Lu2,Ld2,Ru2,Rd2;
    private:
        QCDOrder qcdorder;
        vector<int> partons {-5,-4,-3,-2,-1,1,2,3,4,5,21};
        map<int,double> SigRcoef;
        std::string pdfname;
        LHAPDF::PDFSet * set;
        LHAPDF::PDFUncertainty xuErr;
        bool is_var;
        int ivar=0;
        size_t nmem;
        vector<LHAPDF::PDF*> pdfs ;

        double d_nucleon, d_lepton;
        double M_lepton = -1.;
        double CP_factor = std::numeric_limits<double>::max();

        //double SigR_Nu_LO(double, double, vector<LHAPDF::PDFUncertainty>, vector<vector<double>>, int);
        double SigR_Nu_LO(double, double, map<int,LHAPDF::PDFUncertainty>, map<pair<int,int>,double>, int);
        double SigR_Nu_LO(double, double, map<int,double>);
        double SigR_Nu_LO_NC(double, double, map<int,double>);
        double SigR_Nu_LO_NC(double x, double y, map<int,LHAPDF::PDFUncertainty> dis, std::map<std::pair<int,int>,double> cov_m, int c);
        double Evaluate(double, double, double, int);
        double Evaluate(double, double, double);
        double EvaluateVar(double, double, double, int);

        map<int,double> PDFExtract(double, double);

        double ENU = -1;
        bool ienu = false;
        bool INT_TYPE;
        double Mw2,Mz2, M_iso, GF2;
        double M_boson2;
        double error_band;
        map<int,int> parton_num;

        double Xi(double, double);
        double R(double, double);

        template<class T,double (T::*f)(double,double),int n,int m>
        double HGeneric(double, double);
    public:
        PhysConst* pc;

        // us been very bad people
        double Y,X,Q2;
        // =================
        LHAXS();
        LHAXS(string);
        unsigned int GetNumVar() const {
          return pdfs.size();
        }


        double KernelXS_dsdy(double);
        double KernelXS_dsdyVar(double x);
        double dsdy(double);
        double dsdyVar(double y);

        double KernelXS(double*);
        double KernelXS(double*,int);
        double KernelXSVar(double*);
        double KernelXS_TMC(double*);
        template<double (LHAXS::*f)(double*)>
        double VegasIntegratorXS();

        void Set_M_Lepton(double);
        void Set_CP_factor(double);
        void Set_InteractionType(Current);
        void Set_Neutrino_Energy(double);
        void Set_QCDOrder(QCDOrder);
        void Set_Variant(int);

        virtual double F1(double, double);
        virtual double F2(double, double);
        virtual double F3(double, double);
        virtual double xF3(double, double);
        virtual double F4(double, double);
        virtual double F5(double, double);

        double H1(double, double);
        double H2(double, double);
        double H3(double, double);
        double H4(double, double);
        double H5(double, double);

        double G2(double, double);

        double F1_TMC(double, double);
        double F2_TMC(double, double);
        double F3_TMC(double, double);
        double xF3_TMC(double, double);
        double F4_TMC(double, double);
        double F5_TMC(double, double);
        double FL_TMC(double, double);

        //virtual double F1(map<int, double>);
        virtual double F2(map<int, double>&) const;
        virtual double xF3(map<int, double>&) const;
        double SigRed_Evaluate(double, double, double);
        double SigRed_TMC(double, double, double);

        double total();
        double totalVar();

        ~LHAXS(){
            delete pc;
        }
};

//==================================================================================
// VEGAS
//==================================================================================


// real vegas

template<double (LHAXS::*f)(double*)>
double LHAXS::VegasIntegratorXS(){
  if (M_lepton < 0.){
      cerr << "Check lepton mass!" << std::endl;
      exit(-1);
  }

  if (!(CP_factor == 1. || CP_factor == -1)){
      cerr << "Check CP factor!" << std::endl;
      exit(-1);
  }
  if (ENU == -1){
      cerr << "Neutrino energy not set!" << std::endl;
      exit(-1);
  }
  double res,err;
  const unsigned long dim = 2; int calls = 50000;
  double xl[dim] = {  0.0 , 0.0 };
  double xu[dim] = {  1.0 , 1.0 };

  d_nucleon = M_iso / (2. * ENU);
  d_lepton  = SQ(M_lepton) / (2. * M_iso * ENU);

  gsl_rng_env_setup ();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc (T);

  gsl_monte_function F = { &KernelHelper<LHAXS,f>, dim, this};
  gsl_monte_vegas_state *s_vegas = gsl_monte_vegas_alloc (dim);

  // training
  std::cout << "s_vegas: " << s_vegas << std::endl;
  std::cout << &xl << " " << &xu << " "  << dim << " "  << calls << " "  << &r << std::endl;
  gsl_monte_vegas_integrate (&F, xl, xu, dim, 10000, r, s_vegas,
                              &res, &err);
  std::cout << "stop" << std::endl;
  do
  {
  std::cout << "Here: "<<gsl_monte_vegas_chisq (s_vegas) << std::endl;
  gsl_monte_vegas_integrate (&F, xl, xu, dim, calls, r, s_vegas,
                              &res, &err);
//  std::cout << "ChiSq: "<<gsl_monte_vegas_chisq (s_vegas) << std::endl;
  }
  while (fabs (gsl_monte_vegas_chisq (s_vegas) - 1.0) > 0.5 );

  gsl_monte_vegas_free (s_vegas);
  gsl_rng_free (r);
  //std::cout << "Result: " << res << std::endl;
  return res;
}

#endif
