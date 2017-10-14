#include "lhapdf_cross_section.h"

#define SAVE_PATH "./"

int main(int argc, char* argv[]){

  if(argc != 2){
    cerr << "Argument number not valid! Given: " <<  argc << endl;
    cerr << "Usage: pdfname " << endl;
    return 1;
  }

  PhysConst * pc = new PhysConst();

  string pdfname = (string) argv[1];
  LHAXS xs_obj(pdfname);

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

  for (Current IT : {CC,NC}) {
    xs_obj.Set_InteractionType(IT);
    for (NeutrinoType neutype : {neutrino,antineutrino}){
      xs_obj.Set_CP_factor(CP_factor[neutype]);
      for (PDFVar pdfvar : {central}){
        std::string filename_dsdxdy = static_cast<std::string>(SAVE_PATH) + "dsdxdy-"+NeutrinoTypeLabel[neutype]+"-N-"+IntTypeLabel[IT]+"-"+pdfname+"_"+PDFVarLabel[pdfvar]+".dat";

        ofstream outputfile_dsdxdy(filename_dsdxdy.c_str());

        for (double logenu=0;logenu<=7.;logenu+=0.05){
          double enu = pow(10, logenu);
          xs_obj.Set_Neutrino_Energy(enu*pc->GeV);
          for (double logx=-5.;logx<0.;logx+=0.025){
            double x = pow(10, logx);
            for (double logy=-5.;logy<0.;logy+=0.025){
                double y = pow(10., logy);
                double zz[2];
                zz[0] = log(x);
                zz[1] = log(y);

                double dsigdxdy = xs_obj.KernelXS(zz,PDFVarIndex[pdfvar])/cm2;
                outputfile_dsdxdy << enu << "\t"<< x <<  "\t" << y << "\t" << dsigdxdy << std::endl;
            }
          }
        }
        outputfile_dsdxdy.close();
      }
    }
  }

  return 0;
}
