#include "lhapdf_cross_section.h"

//#define SINGLE_ENERGY_TEST
#define SAVE_PATH "/home/carguelles/NuXSSplMkr/data/newxs/"

int main(int argc, char* argv[]){

  if(argc != 3 or argc != 5){
    cerr << "Argument number not valid! Given: " <<  argc << endl;
    cerr << "Usage: pdfname enu" << endl;
    cerr << "Usage: pdfname enu x y" << endl;
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
  //xs_obj.Set_M_Lepton(1.776*xs_obj.pc->GeV);
  xs_obj.Set_M_Lepton(0.105*xs_obj.pc->GeV);

  double cm2 = SQ(pc->cm);
  double enu = atof(argv[2]);
  xs_obj.Set_Neutrino_Energy(enu*pc->GeV);

  double x,y;
  if(argc==5){
    x = atof(argv[3]);
    y = atof(argv[4]);
  }

  std::cout << "Energy: " << enu/pc->GeV << " GeV" << std::endl;
  for (Current IT : {CC,NC}) {
    xs_obj.Set_InteractionType(IT);
    std::cout << "Interaction: " << IntTypeLabel[IT] << std::endl;
    for (NeutrinoType neutype : {neutrino,antineutrino}){
      xs_obj.Set_CP_factor(CP_factor[neutype]);
      std::cout << "NeutrinoType: " << NeutrinoTypeLabel[neutype] << std::endl;
      for (PDFVar pdfvar : {central}){
        std::cout << "PDFStyle: " << PDFVarLabel[pdfvar] << std::endl;
      //for (PDFVar pdfvar : {central,minus,plus}){
        xs_obj.Set_Variant(PDFVarIndex[pdfvar]);
        double sigma = xs_obj.totalVar()/cm2;
        std::cout << "sigma: " << sigma << " cm^2" << std::endl;
        if(argc==5){
          double dsigdy = xs_obj.dsdyVar(y)/cm2;
          std::cout << "dsigmady: " << dsigdy << " cm^2" << std::endl;
        }
      }
    }
  }

  return 0;
}
