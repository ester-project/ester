#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include <cmath>
#include <string.h>
#include "matrix.h"
#include "constants.h"
#include "physics.h"

#include <iostream>
#include <string>
#include <cstring> // For std::strncpy
#include <unistd.h> // for getcwd()

using namespace std;

// Function to get the absolute path
std::string absolute_path(const std::string &relative_path) {
    // Buffer to hold the current working directory
    char cwd[256];

    // Get the current working directory
    if (getcwd(cwd, sizeof(cwd)) == NULL) {
        std::cerr << "Error getting current working directory." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string cwd_str(cwd);

    // Find the position of the directory "Ester" (case-insensitive, so it accounts for ESTER, ester variations)
    std::string base_dir = "Ester";
    std::string::size_type pos = cwd_str.find(base_dir);

    if (pos != std::string::npos) {
        // Include "Ester" and truncate any trailing directories
        std::string path_up_to_base = cwd_str.substr(0, pos + base_dir.length());

        // Construct the full path
        std::string full_path = path_up_to_base + "/" + relative_path;

        return full_path;
    } else {
        std::cerr << "Directory 'Ester' (case-insensitive) not found in the current working directory." << std::endl;
        exit(EXIT_FAILURE); // Or handle the error as needed
    }
}

extern"C" {

    struct PathData {
        char full_path_cstr[256];
        int length;
    };

    void opa_opmesa_(double *, double *,
            double *, double *, double *, double *, double *, double *, double *, char *, size_t,PathData *path_data); 

}

int opa_opmesa(const matrix& X, double Z, const matrix& T, const matrix& rho,
		opa_struct& opa, const double Xsol, const double Ysol, const double Zsol) {
		
    AbundanceMap& abundance_map = global_abundance_map;
    
    int abund_name_length = abundance_map.mixture_name.length(); 
    
    std:string full_path = abundance_map.ester_home;
    
    //std::cout << "abundance_map.ester_home: " << abundance_map.ester_home << std::endl;
    
    int full_path_cstr_length = full_path.length();
    
    
    //std::string relative_path = "tables/op_mono/";
    //std::string full_path = absolute_path(relative_path);

    //int full_path_length = full_path.length(); 
		
    int error = 0;
    static bool init = false;

	opa.k.dim(T.nrows(), T.ncols());
    opa.xi.dim(T.nrows(), T.ncols());
    opa.dlnxi_lnrho.dim(T.nrows(), T.ncols());
    opa.dlnxi_lnT.dim(T.nrows(), T.ncols());

   
    for (int i=0; i<X.nrows(); i++) {
        for (int j=0; j<X.ncols(); j++) {

            //double x[6];
            double x[7];
            double t, ro, kap, dkapt, dkapro, dkapx;
            double dlnkT;
            double abund[17]; 
            double a_weights[17];
            
           char comp_name_cstr[abund_name_length+1]; // Ensure this matches or exceeds the length expected by Fortran
           std::strncpy(comp_name_cstr, abundance_map.mixture_name.c_str(), abund_name_length);
           comp_name_cstr[abund_name_length] = '\0'; // null termination
              
            
           // same as above, just repeadted for full path
           //char full_path_cstr[full_path_cstr_length+1];
           //std::strncpy(full_path_cstr, abundance_map.ester_home.c_str(), full_path_cstr_length);
           //full_path_cstr[full_path_cstr_length] = '\0'; // null termination    
           
           const size_t full_path_cstr_length = 256;  // Fixed length for the path
	   char full_path_cstr[full_path_cstr_length+1];  // +1 for the null terminator
	   std::strncpy(full_path_cstr, abundance_map.ester_home.c_str(), full_path_cstr_length);
           full_path_cstr[full_path_cstr_length] = '\0'; // Null-terminate the string
       
           //std::cout << "C++ full path: " << full_path_cstr << " , Len: " << full_path_cstr_length << endl;          
           
           PathData path_data;
           //path_data.full_path_cstr = full_path_cstr;
           //path_data.length = strlen(path_data.full_path_cstr);	
   	   std::strncpy(path_data.full_path_cstr, full_path_cstr, full_path_cstr_length);  // Copy string into structure
           path_data.length = std::strlen(path_data.full_path_cstr);  // Set the length


           //std::cout << "C++ path_data.full_path: " << path_data.full_path_cstr << " , path_data.len: " << path_data.length << endl;          
               	   
            //int full_path_length = full_path.length();
            //char full_path_cstr[full_path_length + 1];
            //std::strncpy(full_path_cstr, full_path.c_str(), full_path_length);
            //full_path_cstr[full_path_length] = '\0'; // Ensure null termination



            x[0] = X(i, j);
            x[1] = 1.0 - X(i, j) - Z; 
            x[2] = Z; 
            x[3] = Xsol;
            x[4] = Ysol;
            x[5] = Zsol;
            x[6] = abundance_map.M_init;

            
            t = T(i, j);
            ro = rho(i, j);
           
    	    // H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Cr, Mn, Fe, and Ni. Elements considered currently for op mono table integration. 
            
            abund[0] = abundance_map.comp_abund["H"];
            abund[1] = abundance_map.comp_abund["He"];
            abund[2] = abundance_map.comp_abund["C"];
            abund[3] = abundance_map.comp_abund["N"];
            abund[4] = abundance_map.comp_abund["O"];
            abund[5] = abundance_map.comp_abund["Ne"];
            abund[6] = abundance_map.comp_abund["Na"];
            abund[7] = abundance_map.comp_abund["Mg"];
            abund[8] = abundance_map.comp_abund["Al"];
            abund[9] = abundance_map.comp_abund["Si"];
            abund[10] = abundance_map.comp_abund["S"];
            abund[11] = abundance_map.comp_abund["Ar"];
            abund[12] = abundance_map.comp_abund["Ca"];
            abund[13] = abundance_map.comp_abund["Cr"];
            abund[14] = abundance_map.comp_abund["Mn"];
            abund[15] = abundance_map.comp_abund["Fe"];
            abund[16] = abundance_map.comp_abund["Ni"];         
            
            a_weights[0] = abundance_map.A_weights["H"];
            a_weights[1] = abundance_map.A_weights["He"];
            a_weights[2] = abundance_map.A_weights["C"];
            a_weights[3] = abundance_map.A_weights["N"];
            a_weights[4] = abundance_map.A_weights["O"];
            a_weights[5] = abundance_map.A_weights["Ne"];
            a_weights[6] = abundance_map.A_weights["Na"];
            a_weights[7] = abundance_map.A_weights["Mg"];
            a_weights[8] = abundance_map.A_weights["Al"];
            a_weights[9] = abundance_map.A_weights["Si"];
            a_weights[10] = abundance_map.A_weights["S"];
            a_weights[11] = abundance_map.A_weights["Ar"];
            a_weights[12] = abundance_map.A_weights["Ca"];
            a_weights[13] = abundance_map.A_weights["Cr"];
            a_weights[14] = abundance_map.A_weights["Mn"];
            a_weights[15] = abundance_map.A_weights["Fe"];
            a_weights[16] = abundance_map.A_weights["Ni"];    
            
            
            //opa_opmesa_(x, &t, &ro, &kap, &dkapt, &dkapro, &dkapx, abund, a_weights, comp_name_cstr,abund_name_length,full_path_cstr,full_path_length);
            
            
            //opa_opmesa_(x, &t, &ro, &kap, &dkapt, &dkapro, &dkapx, abund, a_weights, comp_name_cstr,abund_name_length);
                        
                        //opa_opmesa_(x, &t, &ro, &kap, &dkapt, &dkapro, &dkapx, abund, a_weights, comp_name_cstr,abund_name_length,full_path_cstr,full_path_cstr_length);
            
            
            opa_opmesa_(x, &t, &ro, &kap, &dkapt, &dkapro, &dkapx, abund, a_weights, comp_name_cstr, abund_name_length, &path_data);

            
            //opa_opmesa_(x, &t, &ro, &kap, &dkapt, &dkapro, &dkapx, abund, a_weights, comp_name_cstr,abund_name_length,full_path_cstr,full_path_cstr_length);            
            
            
            
            //opa_opmesa_(x, &t, &ro, &kap, &dkapt, &dkapro, &dkapx, abund, a_weights, comp_name_cstr,abund_name_length,full_path_cstr);
            
            opa.k(i, j)=pow(10,kap);
 	        dlnkT-=3*dkapro;
	        opa.xi(i, j)=16*SIG_SB*pow(t,3)/(3*opa.k(i, j)*ro);
	        opa.dlnxi_lnrho(i ,j)=(-1-dkapro);
            opa.dlnxi_lnT(i, j)=(3-dkapt);
		

        }

    }
    

    //printf("OP   Centre:  T = %e rho = %e kap = %e, dlnxi_lnT = %e, dlnxi_lnrho = %e \n", T(1), rho(1), opa.k(1), opa.dlnxi_lnT(1), opa.dlnxi_lnrho(1));
    //printf("OP   Surface: T = %e rho = %e kap = %e, dlnxi_lnT = %e, dlnxi_lnrho = %e \n", T(-1), rho(-1), opa.k(-1), opa.dlnxi_lnT(-1), opa.dlnxi_lnrho(-1));

    return error;

}
