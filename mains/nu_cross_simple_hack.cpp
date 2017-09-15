#include "lhapdf_cross_section.h"

//#define SINGLE_ENERGY_TEST
#define SAVE_PATH "/home/carguelles/NuXSSplMkr/data/test/"

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

  //enum IntType {CC,NC};
  enum NeutrinoType {neutrino,antineutrino};
  enum PDFVar {central,minus,plus};

  std::map<Current,std::string> IntTypeLabel {{CC,"cc"},{NC,"nc"}};
  std::map<NeutrinoType,double> CP_factor {{neutrino,1.},{antineutrino,-1}};
  std::map<NeutrinoType,std::string> NeutrinoTypeLabel {{neutrino,"numu"},{antineutrino,"numubar"}};
  std::map<PDFVar,int> PDFVarIndex {{central,0},{minus,-1},{plus,1}};
  std::map<PDFVar,std::string> PDFVarLabel {{central,"central"},{minus,"minus"},{plus,"plus"}};

  // muon mass
  xs_obj.Set_M_Lepton(0.105*xs_obj.pc->GeV);

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
    for (PDFVar pdfvar : {central,minus,plus}){
      double dsigdxdy = xs_obj.KernelXS(ZZ,PDFVarIndex[pdfvar])/cm2;
      std::cout << ENU << "\t"<< ZZ[0] <<  "\t" << ZZ[1] << "\t";
      std::cout << PDFVarLabel[pdfvar] << "\t" << dsigdxdy << std::endl;
    }
  }
  exit(1);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
#endif

  for (Current IT : {CC,NC}) {
    xs_obj.Set_InteractionType(IT);
    for (NeutrinoType neutype : {neutrino,antineutrino}){
      xs_obj.Set_CP_factor(CP_factor[neutype]);
      for (PDFVar pdfvar : {central}){
        std::string filename_dsdy = static_cast<std::string>(SAVE_PATH) + "dsde-"+NeutrinoTypeLabel[neutype]+"-N-"+IntTypeLabel[IT]+"-"+pdfname+"_"+PDFVarLabel[pdfvar]+".dat";
        ofstream outputfile_dsdy(filename_dsdy.c_str());
        for (double logenu=0;logenu<=7.;logenu+=0.05){
          double enu = pow(10, logenu);
          xs_obj.Set_Neutrino_Energy(enu*pc->GeV);
          for (double logenu_bis=0;logenu_bis<=7.;logenu_bis+=0.05){
            double enu_bis = pow(10,logenu_bis);
            double y = 1. - enu_bis/enu;
            if(y>1 or y <0)
              outputfile_dsdy << enu << "\t" << enu_bis << "\t"<< 0. << std::endl;
            else{
              double dsigdy = xs_obj.dsdy(y)/cm2;
              outputfile_dsdy << enu << "\t" << enu_bis << "\t"<< dsigdy/enu << std::endl;
            }
          }
        }
        outputfile_dsdy.close();
      }
    }
  }

  return 0;
}
