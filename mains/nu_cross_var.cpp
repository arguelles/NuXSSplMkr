#include "lhapdf_cross_section.h"

//#define SINGLE_ENERGY_TEST
#define SAVE_PATH "/home/carguelles/NuXSSplMkr/data/VAR/"

int main(int argc, char* argv[]){

#ifndef SINGLE_ENERGY_TEST
  if(argc != 2){
    cerr << "Argument number not valid! Given: " <<  argc << endl;
    cerr << "Usage: pdfname " << endl;
    return 1;
  }
#else
  if(argc != 5){
    cerr << "Argument number not valid! Given: " <<  argc << endl;
    cerr << "Usage: pdfname enu x y" << endl;
    return 1;
  }
#endif

  PhysConst * pc = new PhysConst();

  string pdfname = (string) argv[1];
  LHAXS xs_obj(pdfname);

  enum NeutrinoType {neutrino,antineutrino};
  enum PDFVar {central,minus,plus};

  std::map<NeutrinoType,double> CP_factor {{neutrino,1.},{antineutrino,-1}};
  std::map<NeutrinoType,std::string> NeutrinoTypeLabel {{neutrino,"numu"},{antineutrino,"numubar"}};

  // muon mass
  xs_obj.Set_M_Lepton(0.105*xs_obj.pc->GeV);
  // set CC interaction
  xs_obj.Set_InteractionType(CC);

  double cm2 = SQ(pc->cm);

#ifdef SINGLE_ENERGY_TEST
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//158.489 0.891251  0.1
  double ENU = atof(argv[2]);
  xs_obj.Set_Neutrino_Energy(ENU*pc->GeV);
  double ZZ[2];
  ZZ[0] = atof(argv[3]);
  ZZ[1] = atof(argv[4]);

  for (NeutrinoType neutype : {neutrino,antineutrino}){
    std::cout << "==== " << NeutrinoTypeLabel[neutype] << " ====" << std::endl;
    xs_obj.Set_CP_factor(CP_factor[neutype]);
    std::cout << ENU << "\t"<< ZZ[0] <<  "\t" << ZZ[1] << "\t";
    for (unsigned int var = 0; var < xs_obj.GetNumVar(); var++){
    //for (unsigned int var = 0; var < xs_obj.GetNumVar(); var++){
      double dsigdxdy = xs_obj.KernelXSVar(ZZ)/cm2;
      std::cout << dsigdxdy << "\t";
    }
    std::cout << std::endl;
  }
  exit(1);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
#endif

  for (NeutrinoType neutype : {neutrino,antineutrino}){
      xs_obj.Set_CP_factor(CP_factor[neutype]);
      // neutrino = 1., antineutrino = -1
      std::string filename_dsdxdy = static_cast<std::string>(SAVE_PATH) + "dsdxdy-"+NeutrinoTypeLabel[neutype]+"-N-cc-"+pdfname+"_variations.dat";
      std::string filename_dsdy = static_cast<std::string>(SAVE_PATH) + "dsdy-"+NeutrinoTypeLabel[neutype]+"-N-cc-"+pdfname+"_variations.dat";

      ofstream outputfile_dsdxdy(filename_dsdxdy.c_str());
      ofstream outputfile_dsdy(filename_dsdy.c_str());

      for (double logenu=2;logenu<=6.;logenu+=0.05){
        double enu = pow(10, logenu);
        xs_obj.Set_Neutrino_Energy(enu*pc->GeV);
        for (double logx=-5.;logx<0.;logx+=0.025){
          double x = pow(10, logx);
          for (double logy=-5.;logy<0.;logy+=0.025){
              double y = pow(10, logy);
              double zz[2];
              zz[0] = x;
              zz[1] = y;
              outputfile_dsdxdy << enu << "\t"<< x <<  "\t" << y << "\t";
              for (unsigned int var = 0; var < xs_obj.GetNumVar(); var++){
                 double dsigdxdy = xs_obj.KernelXSVar(zz)/cm2;
                 outputfile_dsdxdy << dsigdxdy << "\t";
              }
              outputfile_dsdxdy << std::endl;
          }
          //double dsigdy = xs_obj.dsdy(x,PDFVarIndex[pdfvar])/cm2;
          double dsigdy = 0.;
          //outputfile_dsdy << enu << "\t"<< x << "\t" << dsigdy << std::endl;
        }
      }

      outputfile_dsdxdy.close();
      outputfile_dsdy.close();
  }

  return 0;
}
