#include "lhapdf_cross_section.h"

//#define SINGLE_ENERGY_TEST
#define SAVE_PATH "/home/carguelles/NuXSSplMkr/data/newxs/"

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

         std::vector<double> e_range {1.00000000e+03,   1.08436597e+03,   1.17584955e+03,
         1.27505124e+03,   1.38262217e+03,   1.49926843e+03,
         1.62575567e+03,   1.76291412e+03,   1.91164408e+03,
         2.07292178e+03,   2.24780583e+03,   2.43744415e+03,
         2.64308149e+03,   2.86606762e+03,   3.10786619e+03,
         3.37006433e+03,   3.65438307e+03,   3.96268864e+03,
         4.29700470e+03,   4.65952567e+03,   5.05263107e+03,
         5.47890118e+03,   5.94113398e+03,   6.44236351e+03,
         6.98587975e+03,   7.57525026e+03,   8.21434358e+03,
         8.90735464e+03,   9.65883224e+03,   1.04737090e+04,
         1.13573336e+04,   1.23155060e+04,   1.33545156e+04,
         1.44811823e+04,   1.57029012e+04,   1.70276917e+04,
         1.84642494e+04,   2.00220037e+04,   2.17111795e+04,
         2.35428641e+04,   2.55290807e+04,   2.76828663e+04,
         3.00183581e+04,   3.25508860e+04,   3.52970730e+04,
         3.82749448e+04,   4.15040476e+04,   4.50055768e+04,
         4.88025158e+04,   5.29197874e+04,   5.73844165e+04,
         6.22257084e+04,   6.74754405e+04,   7.31680714e+04,
         7.93409667e+04,   8.60346442e+04,   9.32930403e+04,
         1.01163798e+05,   1.09698580e+05,   1.18953407e+05,
         1.28989026e+05,   1.39871310e+05,   1.51671689e+05,
         1.64467618e+05,   1.78343088e+05,   1.93389175e+05,
         2.09704640e+05,   2.27396575e+05,   2.46581108e+05,
         2.67384162e+05,   2.89942285e+05,   3.14403547e+05,
         3.40928507e+05,   3.69691271e+05,   4.00880633e+05,
         4.34701316e+05,   4.71375313e+05,   5.11143348e+05,
         5.54266452e+05,   6.01027678e+05,   6.51733960e+05,
         7.06718127e+05,   7.66341087e+05,   8.30994195e+05,
         9.01101825e+05,   9.77124154e+05,   1.05956018e+06,
         1.14895100e+06,   1.24588336e+06,   1.35099352e+06,
         1.46497140e+06,   1.58856513e+06,   1.72258597e+06,
         1.86791360e+06,   2.02550194e+06,   2.19638537e+06,
         2.38168555e+06,   2.58261876e+06,   2.80050389e+06,
         3.03677112e+06,   3.29297126e+06,   3.57078596e+06,
         3.87203878e+06,   4.19870708e+06,   4.55293507e+06,
         4.93704785e+06,   5.35356668e+06,   5.80522552e+06,
         6.29498899e+06,   6.82607183e+06,   7.40196000e+06,
         8.02643352e+06,   8.70359136e+06,   9.43787828e+06,
         1.02341140e+07,   1.10975250e+07,   1.20337784e+07,
         1.30490198e+07,   1.41499130e+07,   1.53436841e+07,
         1.66381689e+07,   1.80418641e+07,   1.95639834e+07,
         2.12145178e+07,   2.30043012e+07,   2.49450814e+07,
         2.70495973e+07,   2.93316628e+07,   3.18062569e+07,
         3.44896226e+07,   3.73993730e+07,   4.05546074e+07,
         4.39760361e+07,   4.76861170e+07,   5.17092024e+07,
         5.60716994e+07,   6.08022426e+07,   6.59318827e+07,
         7.14942899e+07,   7.75259749e+07,   8.40665289e+07,
         9.11588830e+07,   9.88495905e+07,   1.07189132e+08,
         1.16232247e+08,   1.26038293e+08,   1.36671636e+08,
         1.48202071e+08,   1.60705282e+08,   1.74263339e+08,
         1.88965234e+08,   2.04907469e+08,   2.22194686e+08,
         2.40940356e+08,   2.61267523e+08,   2.83309610e+08,
         3.07211300e+08,   3.33129479e+08,   3.61234270e+08,
         3.91710149e+08,   4.24757155e+08,   4.60592204e+08,
         4.99450512e+08,   5.41587138e+08,   5.87278661e+08,
         6.36824994e+08,   6.90551352e+08,   7.48810386e+08,
         8.11984499e+08,   8.80488358e+08,   9.54771611e+08,
         1.03532184e+09,   1.12266777e+09,   1.21738273e+09,
         1.32008840e+09,   1.43145894e+09,   1.55222536e+09,
         1.68318035e+09,   1.82518349e+09,   1.97916687e+09,
         2.14614120e+09,   2.32720248e+09,   2.52353917e+09,
         2.73644000e+09,   2.96730241e+09,   3.21764175e+09,
         3.48910121e+09,   3.78346262e+09,   4.10265811e+09,
         4.44878283e+09,   4.82410870e+09,   5.23109931e+09,
         5.67242607e+09,   6.15098579e+09,   6.66991966e+09,
         7.23263390e+09,   7.84282206e+09,   8.50448934e+09,
         9.22197882e+09,   1.00000000e+10 };

  PhysConst * pc = new PhysConst();

  string pdfname = (string) argv[1];
  LHAXS xs_obj(pdfname);

  //enum IntType {CC,NC};
  enum NeutrinoType {neutrino,antineutrino};
  enum PDFVar {central,minus,plus};

  std::map<Current,std::string> IntTypeLabel {{CC,"cc"},{NC,"nc"}};
  std::map<NeutrinoType,double> CP_factor {{neutrino,1.},{antineutrino,-1}};
  std::map<NeutrinoType,std::string> NeutrinoTypeLabel {{neutrino,"nutau"},{antineutrino,"nutaubar"}};
  std::map<PDFVar,int> PDFVarIndex {{central,0},{minus,-1},{plus,1}};
  std::map<PDFVar,std::string> PDFVarLabel {{central,"central"},{minus,"minus"},{plus,"plus"}};

  // muon mass
  xs_obj.Set_M_Lepton(1.776*xs_obj.pc->GeV);

  double cm2 = SQ(pc->cm);

  for (Current IT : {CC,NC}) {
    xs_obj.Set_InteractionType(IT);
    for (NeutrinoType neutype : {neutrino,antineutrino}){
      xs_obj.Set_CP_factor(CP_factor[neutype]);
      //for (PDFVar pdfvar : {central}){
      for (PDFVar pdfvar : {central,minus,plus}){
        std::cout << "BEGIN sigma-"+NeutrinoTypeLabel[neutype]+"-N-"+IntTypeLabel[IT]+"-"+pdfname+"_"+PDFVarLabel[pdfvar]+".dat" << std::endl;
        xs_obj.Set_Variant(PDFVarIndex[pdfvar]);
        // neutrino = 1., antineutrino = -1
        std::string filename_dsdy = static_cast<std::string>(SAVE_PATH) + "dsdE-"+NeutrinoTypeLabel[neutype]+"-N-"+IntTypeLabel[IT]+"-"+pdfname+"_"+PDFVarLabel[pdfvar]+".dat";
        std::string filename_sigma = static_cast<std::string>(SAVE_PATH) + "sigma-"+NeutrinoTypeLabel[neutype]+"-N-"+IntTypeLabel[IT]+"-"+pdfname+"_"+PDFVarLabel[pdfvar]+".dat";

        ofstream outputfile_dsdy(filename_dsdy.c_str());
        ofstream outputfile_sigma(filename_sigma.c_str());

        //for (double logenu=0;logenu<=7.;logenu+=0.05){
        //  double enu = pow(10, logenu);
        for( double enu : e_range){
          xs_obj.Set_Neutrino_Energy(enu*pc->GeV);
          //for (double logenu2=0;logenu2<=7.;logenu2+=0.05){
          //  double enu2= pow(10, logenu2);
          for ( double enu2 : e_range){
            double y = 1. - enu2/enu;
            double dsigdy=0;
            if (y > 0.)
              dsigdy = xs_obj.dsdyVar(y)/cm2;
            outputfile_dsdy << dsigdy/enu << "\t";
          }
          outputfile_dsdy << endl;
          double sigma = xs_obj.totalVar();
          outputfile_sigma << enu << "\t"<< sigma/cm2 << std::endl;
        }

        outputfile_dsdy.close();
        outputfile_sigma.close();
        std::cout << "END sigma-"+NeutrinoTypeLabel[neutype]+"-N-"+IntTypeLabel[IT]+"-"+pdfname+"_"+PDFVarLabel[pdfvar]+".dat" << std::endl;
      }
    }
  }

  return 0;
}
