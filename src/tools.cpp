#include "tools.h"

std::string toString(double value) {
   std::stringstream ss;
   ss << value;
   return ss.str();
}

std::string toString(int value) {
   std::stringstream ss;
   ss << value;
   return ss.str();
}


bool fexists(string filepath)
{
  ifstream ifile(filepath.c_str());
  return ifile;
};

Table quickread(string filepath){
    // create and open file stream
    ifstream infile(filepath.c_str());

    if(!infile){
        cerr << "Error: file could not be opened. Filepath " << filepath.c_str()<<  endl;
        exit(0);
    }
    #ifdef quickread_DEBUG
    int x,y;
    #endif
    Table table;

    string line;
    while(getline(infile,line)){
        Row row;
        stringstream linestream(line);

        double data;
        while(linestream >> data){
            row.push_back(data);
        }
        if (!row.empty()){
            #ifdef quickread_DEBUG
            y = row.size();
            #endif
            table.push_back(row);
        }
    }

    #ifdef quickread_DEBUG
    x = table.size();
    cout << "x: " << x << " y: "<< y << endl;
    cout << table[10][0] << endl;
    #endif
    return table;
};

int quickread(string filepath,Table tbl){
    // create and open file stream
    ifstream infile(filepath.c_str());

    if(!infile){
        cerr << "Error: file could not be opened. Filepath " << filepath.c_str()<<  endl;
        exit(0);
    }
    #ifdef quickread_DEBUG
    int x,y;
    #endif

    string line;
    while(getline(infile,line)){
        Row row;
        stringstream linestream(line);

        double data;
        while(linestream >> data){
            row.push_back(data);
        }
        if (!row.empty()){
            #ifdef quickread_DEBUG
            y = row.size();
            #endif
            tbl.push_back(row);
        }
    }

    #ifdef quickread_DEBUG
    x = tbl.size();
    cout << "x: " << x << " y: "<< y << endl;
    cout << tbl[10][0] << endl;
    #endif
    return 0;
};

int quickwrite(string filepath, Table tbl){
    // create and open file stream
    ofstream outfile(filepath.c_str());

    if(!outfile){
        cerr << "Error: file could not be created. Filepath " << filepath.c_str()<<  endl;
        exit(0);
    }

    int tbl_size = tbl.size();

    for (int i=0; i < tbl_size;i++){
        vector<double> line = tbl[i];
        int line_size = line.size();
        for(int j=0; j < line_size; j++){
           outfile << line[j] << " ";
        }
        outfile << endl;
    }

    outfile.close();

    return 0;
};

vector<double> linspace(double Emin,double Emax,int div){
    vector<double> linpoints;
    double step_lin = (Emax - Emin)/double(div);

    double EE = Emin;
    while (EE <= Emax+0.001){
        linpoints.push_back(EE);
        EE = EE + step_lin;
    }
    return linpoints;
};

vector<double> logspace(double Emin,double Emax,int div){
    vector<double> logpoints;
    double Emin_log,Emax_log;
    if (Emin < 1.0e-5 ) {
        Emin_log = 0.0;
    } else {
        Emin_log = log(Emin);
    }
    Emax_log = log(Emax);

    double step_log = (Emax_log - Emin_log)/double(div);

    double EE = Emin_log;
        logpoints.push_back(Emin_log);
    while (EE <= Emax_log+0.001){
        logpoints.push_back(exp(EE));
        EE = EE + step_log;
    }
    return logpoints;
};

//Table intertable(Table xy_data, vector<double> x_array, int j1 = 0, int j2 = 1){
//    Table result;
//    int arraysize = xy_data.size();
//    
//    double xx[arraysize];
//    double yy[arraysize];
//            
//    for (int i=0; i < arraysize;i++){
//        xx[i] = xy_data[i][j1];
//        yy[i] = xy_data[i][j2];
//    }
//            
//    gsl_spline* inter = gsl_spline_alloc(gsl_interp_cspline,arraysize);
//    gsl_interp_accel* inter_accel = gsl_interp_accel_alloc ();
//    gsl_spline_init(inter,xx,yy,arraysize);
//    
//    for(unsigned int i=0; i < x_array.size();i++){
//        Row row;
//        row.push_back(x_array[i]);
//        if (x_array[i] > xx[arraysize-1] or x_array[i] < xx[0]){
//            row.push_back(0.0);
//        } else {
//            row.push_back(gsl_spline_eval(inter,x_array[i],inter_accel));    
//        }
//        result.push_back(row);
//    }
//    
//    gsl_spline_free(inter);
//    gsl_interp_accel_free(inter_accel);
//    
//    return result;
//};

void PrintTable(Table tbl){
    for(int i = 0; i < tbl.size(); i++){
        for(int j = 0; j < tbl[i].size(); j++){
            std::cout << tbl[i][j] <<  " ";
        }
        std::cout << std::endl;
    }
}

double Abs(double x){
    return (x<0) ? -x : x;
}

double FindNearestPointInVector(vector<double> k, double x)
{
    double point = k[0];
    //double min = abs ( x - k[0] );
    double min =  Abs (x - k[0]) ;
    for (int i = 0; i<k.size(); i++){
        double diff = Abs ( x - k[i] );
        if (diff == min){
            continue;
        }
        if (diff < min){
            min = diff;
            point = k[i];
        }
    }

    return point;
}

vector< vector<double> > TableToVector(string path){
    double a,b,c,d,e,f;
    vector<double> va,vb,vc,vd,ve,vf;
    vector< vector<double> > complete;
    std::ifstream file(path.c_str());
    
    // Reading columns into vectors
    while(!file.eof()){
        file >> a >> b >> c >> d >> e >> f;
        va.push_back(a);
        vb.push_back(b);
        vc.push_back(c);
        vc.push_back(d);
        vc.push_back(e);
        vc.push_back(f);
    }

    // Delete last element
    va.pop_back();
    vb.pop_back();
    vc.pop_back();
    vd.pop_back();
    ve.pop_back();
    vf.pop_back();

    // Fill in to the complete vector
    complete.push_back(va);
    complete.push_back(vb);
    complete.push_back(vc);
    complete.push_back(vd);
    complete.push_back(ve);
    complete.push_back(vf);

    return complete;

}

// additional GSL-like tools

//void gsl_matrix_complex_conjugate(gsl_matrix_complex *cm)
//{
//  gsl_complex z;
//  size_t i, j;
//  for (i = 0; i < cm->size1; i++) {
//    for (j = 0; j < cm->size2; j++) {
//      z = gsl_matrix_complex_get(cm, i, j);
//      gsl_matrix_complex_set(cm, i, j, gsl_complex_conjugate(z));
//    }
//  }
//};
//
//void gsl_matrix_complex_print(gsl_matrix_complex* matrix){
//    for(unsigned int i = 0; i < matrix->size1; i++){
//        for(unsigned int j = 0; j < matrix->size2; j++){
//            cout << gsl_matrix_complex_get(matrix,i,j).dat[0] <<
//            "+i" << gsl_matrix_complex_get(matrix,i,j).dat[1] << " ";
//        }
//        cout << endl;
//    }
//};
//
//void gsl_matrix_complex_change_basis_UMUC(gsl_matrix_complex* U, gsl_matrix_complex* M){
//    int numneu = U->size1;
//    gsl_matrix_complex *U1 = gsl_matrix_complex_alloc(numneu,numneu);
//    gsl_matrix_complex *U2 = gsl_matrix_complex_alloc(numneu,numneu);
//    gsl_matrix_complex_memcpy(U1,U);
//    gsl_matrix_complex_memcpy(U2,U);
//    
//    gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(numneu,numneu);
//    
//    // doing : U M U^dagger
//    
//    gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,
//                   gsl_complex_rect(1.0,0.0),M,
//                   U1,gsl_complex_rect(0.0,0.0),T1);
//    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,
//                   gsl_complex_rect(1.0,0.0),U2,
//                   T1,gsl_complex_rect(0.0,0.0),M);
//    // now H_current is in the interaction basis of the mass basis
//    
//    gsl_matrix_complex_free(U1);
//    gsl_matrix_complex_free(U2);
//    gsl_matrix_complex_free(T1);
//};
//
//void gsl_matrix_complex_change_basis_UCMU(gsl_matrix_complex* U, gsl_matrix_complex* M){
//    int numneu = U->size1;
//    gsl_matrix_complex *U1 = gsl_matrix_complex_alloc(numneu,numneu);
//    gsl_matrix_complex *U2 = gsl_matrix_complex_alloc(numneu,numneu);
//    gsl_matrix_complex_memcpy(U1,U);
//    gsl_matrix_complex_memcpy(U2,U);
//    
//    gsl_matrix_complex *T1 = gsl_matrix_complex_alloc(numneu,numneu);
//    
//    // doing : U M U^dagger
//    
//    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,
//                   gsl_complex_rect(1.0,0.0),M,
//                   U1,gsl_complex_rect(0.0,0.0),T1);
//    gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,
//                   gsl_complex_rect(1.0,0.0),U2,
//                   T1,gsl_complex_rect(0.0,0.0),M);
//    // now H_current is in the interaction basis of the mass basis
//    
//    gsl_matrix_complex_free(U1);
//    gsl_matrix_complex_free(U2);
//    gsl_matrix_complex_free(T1);
//};

double LinInter(double x,double xM, double xP,double yM,double yP){
    return yM + (yP-yM)*(x-xM)/(xP-xM);
};
