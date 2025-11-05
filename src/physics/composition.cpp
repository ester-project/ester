#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "physics.h"

#include <iostream>
#include <unistd.h>
#include <fstream>

#include <vector> 
#include <sstream>

//#include <filesystem> // For std::filesystem
#include <unistd.h> // For getcwd function on Unix-like systems
#include <cstring>
#include <dirent.h> // For directory traversal on Unix-like systems

#include <string>
#include <map>
#include <cmath>

//#include <boost/stacktrace.hpp>
#include "star.h" // adding for abundance name 

#include <iomanip>

std::string GetCurrentWorkingDirectory() {
    const size_t bufferSize = 1024;
    char buffer[bufferSize];
    if (getcwd(buffer, bufferSize) != nullptr) {
        return std::string(buffer);
    }
    return std::string(""); 
}

//------
// Define the logging function and macro in main.cpp
/*inline void logCaller(const std::string& functionName, const std::string& fileName, int lineNumber) {
    std::cout << "Called from function: " << functionName 
              << " in file: " << fileName 
              << " at line: " << lineNumber << std::endl;
	// Print the stack trace
    std::cout << boost::stacktrace::stacktrace() << std::endl;	

	std::string currentPath = GetCurrentWorkingDirectory();

	std::cout << "CWD stack trace call:" << currentPath << std::endl;
}

#define LOG_CALLER() logCaller(__func__, __FILE__, __LINE__)
*/

std::string cutOffPath(const std::string& path, const std::string& delimiter) {
    size_t pos = path.rfind(delimiter); // Find the last occurrence of the delimiter
    if (pos != std::string::npos) {
        return path.substr(0, pos + delimiter.size());
    }
    return path; // If the delimiter is not found, return the original path
}


std::string readPathFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return "";
    }

    std::string path;
    if (std::getline(file, path)) {
        return path;
    } else {
        std::cerr << "Error reading path from file: " << filename << std::endl;
        return "";
    }
}


// creating a structure to fill up element data --> necessary? Only used once. 

struct ElementData {
    int Z;
    std::string El;
    int A;
    float percentage;
    float N;
    float mass_excess;
};

using namespace std;

// defined in physics.h holds abundances, atomic weights and name of abundance file. 
//AbundanceMap global_abundance_map;

//FileMeta File_meta_data;

// Use this if you are not updating the composition on every step. Otherwise
// Call save the results of parse_composition_data and call only update_initial_composition.
double_map initial_composition(double X, double Z) {
    CompositionData initial_comp = parse_composition_data();
    return update_initial_composition(initial_comp, X, Z);
}

CompositionData parse_composition_data() {

	// The environment variable ESTER should be set to the ESTER directory
    // in .bashrc or .cshrc

    // we read the environment variable ESTER
    std::string esterDirectory=getenv("ESTER");
    global_abundance_map.ester_home = esterDirectory;
    
    /** Read data from lodders03_w_mass_excess.txt and populate the map
    this can be done outside of the initial composition function I believe.
    the initial composition should only be done once, however in the evolution
    branch comp will be updated
    **/

    // load lodders03_data_w_mass_excess file 
    std::string lodders03_comp_name = "lodders03_data_w_mass_excess.txt";
    ifstream file_mass_excess(esterDirectory+"/Solar_compositions/"+lodders03_comp_name);
	if (!file_mass_excess.is_open()) {
		// check that file exists
        printf("-------------------- \n ");
        //printf("ESTER root directory is not $HOME/Documents/Ester\n");
        //printf("please make a change in src/physics/composition.cpp l.120\n");
        cerr << "Error opening file: " << esterDirectory+"/Solar_compositions/"+lodders03_comp_name << endl; 
        printf("-------------------- \n ");
    }

	std::map<std::string, std::vector<ElementData>> comp_inp;// this is the variable which is declared here will be filled up with isotope info
	std::string line_mass_excess;
	int count_mass_excess=0;
	while (std::getline(file_mass_excess, line_mass_excess)) {
		if (count_mass_excess ==0) {
			count_mass_excess+=1;
			continue;
		}
		std::istringstream iss(line_mass_excess);
		ElementData data;
		iss >> data.Z >> data.El >> data.A >> data.percentage >> data.N >> data.mass_excess;

		comp_inp[data.El].push_back(data);
		count_mass_excess+=1;

	}

    std::map<std::string, double> comp_out; // converted abundances stored here to be normalised and used for comp
	std::map<std::string, double> comp_A; // atomic weights which are required for -opa mono mode. 
    string line_inp;

    //AbundanceMap& abundance_map = global_abundance_map;
    std::string abund_inp_comp_name = global_abundance_map.mixture_name; 
    
    ifstream file_abund_inp(esterDirectory+"/Solar_compositions/"+abund_inp_comp_name+"_ESTER_abund_input.txt");
	if (!file_abund_inp.is_open()) {
        cerr << "Error opening file: " << esterDirectory+"/Solar_compositions/"+abund_inp_comp_name+"_ESTER_abund_input.txt"<< endl; 
    	}

    std::map<std::string, float> comp_abund_test; 

	bool in_header = false;

	// Read data from the input abundance file and use each element as "element_to_search"
	while (std::getline(file_abund_inp, line_inp)) {

        // Check for header markers
        if (line_inp.find("### header ###") != std::string::npos) {
            in_header = true;
            continue;
        }

        if (line_inp.find("### end of header ###") != std::string::npos) {
            in_header = false;
            continue;
        }

        // Skip lines if still in header
        if (in_header) {
            continue;
        }

	    // Find the position of any comment character
        size_t comment_pos = line_inp.find_first_of("#!/");
        
        // If found a possible comment start
        if (comment_pos != std::string::npos) {
			// checking for // 
            if (line_inp[comment_pos] == '/' && comment_pos + 1 < line_inp.size() && line_inp[comment_pos + 1] == '/') {
                comment_pos++;  // Skip the second slash bc comment 
            }
            // Trim the string to exclude the comment
            line_inp = line_inp.substr(0, comment_pos);
        }

        
        line_inp.erase(line_inp.find_last_not_of(" \t\n\r\f\v") + 1); // Trim whitespace for clean parsing

        std::string element_to_search;
        float abundance;
        
        // Separate the string and float
        std::stringstream ss(line_inp);
        std::getline(ss, element_to_search, ','); // Extract string before comma
        
        // Trim whitespace around the element name
        element_to_search.erase(0, element_to_search.find_first_not_of(" \t"));
        element_to_search.erase(element_to_search.find_last_not_of(" \t") + 1);

        ss >> abundance; // Extract float after comma

        comp_abund_test[element_to_search]=abundance;

		if (comp_inp.find(element_to_search) != comp_inp.end()) {
			//std::cout << "Data for element " << element_to_search << ":" << std::endl;
			for (const auto& element_data : comp_inp[element_to_search]) {
				// before it was adding element_data.Z, I changed it to A bc isotope, right?
				std::string key = element_to_search + std::to_string(element_data.A);

				// CALCULATION OF ABUNDANCE IN MASS FORM 

				int Z = element_data.Z;  // Number of protons (for example, carbon)
				int N = element_data.A-element_data.Z;  // Number of neutrons (for example, carbon)
				float mass_excess_mev_per_c2 = element_data.mass_excess;  // Mass excess in MeV/c^2 (for example, carbon)

				// Constants
				const float c = 2.998e10;  // Speed of light in vacuum in cm/s

				// Calculate conversion factor from MeV/c^2 to amu in cgs units
				float conversion_factor = 1.602176634e-6 / (c * c) / 1.66053906660e-24;

				// Calculate isotope mass in AMU 
				float isotopic_mass_amu = Z + N + mass_excess_mev_per_c2 * conversion_factor;
				
				//Convert abundance to by mass

				float abundance_mass = std::pow(10,abundance-12) * isotopic_mass_amu * element_data.percentage/100;

				comp_out[key] = abundance_mass;

				comp_A[key] = isotopic_mass_amu * element_data.percentage/100;
			}
		} else {
		}
	}

    std::vector<std::string> opa_list = {"H","He","C", "N", "O", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "Ca", "Ti", "Cr", "Mn", "Fe", "Ni"};
	
	global_abundance_map.comp_abund = comp_abund_test;

	double_map combined_A;
    for (const auto& kv : comp_A) {
        std::string key = kv.first;
        double value = kv.second;

        // Find the index where the numeric part starts
        std::size_t pos = key.find_first_of("0123456789");
        if (pos != std::string::npos) {
            std::string id_part = key.substr(0, pos);  // Extract the ID part (e.g., "Cd" from "Cd111")
            
            // Check if the ID part matches any element in opa_list
            for (const auto& prefix : opa_list) {
                if (id_part == prefix) {
                    // If there is a match, add the value to the corresponding key in the result map
                    if (combined_A.find(prefix) == combined_A.end()) {
                        combined_A[prefix] = 0.0;
                    }
                    combined_A[prefix] += value;
                                     
                }
            }
        }
    }

    double z_sum = 0.0;
	double y_sum = 0.0;
	double x_sum = 0.0;
    std::map<std::string, double> comp_out_x;
    std::map<std::string, double> comp_out_y;
    std::map<std::string, double> comp_out_z;

	for (const auto& entry : comp_out) {
        const std::string& key = entry.first;
        double value = entry.second;

        if (key.find("H1") != std::string::npos || key.find("H2") != std::string::npos) {
            comp_out_x[key] = value;
            x_sum += value;
        } else if (key.find("He3") != std::string::npos || key.find("He4") != std::string::npos) {
            comp_out_y[key] = value;
            y_sum += value;
        } else {
            comp_out_z[key] = value;
            z_sum += value;
        }
    }
    double xyz_sum = x_sum + y_sum + z_sum;
    
    // Divide all values of comp_out_x by x_sum
    for (auto& entry : comp_out_x) {
        entry.second /= x_sum;
    }
    // Divide all values of comp_out_y by y_sum
    for (auto& entry : comp_out_y) {
        entry.second /= y_sum;
    }
    // Divide all values of comp_out_z by z_sum
    for (auto& entry : comp_out_z) {
        entry.second /= z_sum;
    }
    

    global_abundance_map.A_weights = combined_A;

    std::map<std::string, double> comp_test_xyz = comp_out; // "" by xyz_sum 
    double tot_test_xyz = 0.0;
	double tot_test_xyz_x_only = 0.0;
	double tot_test_xyz_y_only = 0.0;
	double tot_test_xyz_z_only = 0.0;

    for (auto& entry : comp_test_xyz) {
        entry.second /= xyz_sum;
		tot_test_xyz += entry.second;
	
	    const std::string& key = entry.first;
        double value = entry.second;

        if (key.find("H1") != std::string::npos || key.find("H2") != std::string::npos) {
            //comp_out_x[key] = value;
            tot_test_xyz_x_only += value;
        } else if (key.find("He3") != std::string::npos || key.find("He4") != std::string::npos) {
            //comp_out_y[key] = value;
            tot_test_xyz_y_only += value;
        } else {
            //comp_out_z[key] = value;
            tot_test_xyz_z_only += value;
        }
	}

    for (auto& entry : comp_test_xyz) {

		const std::string& key = entry.first;
        double value = entry.second;

        if (key.find("H1") != std::string::npos || key.find("H2") != std::string::npos) {
			entry.second /= tot_test_xyz_x_only; 
        } else if (key.find("He3") != std::string::npos || key.find("He4") != std::string::npos) {
			entry.second /= tot_test_xyz_y_only;
        } else {
			entry.second /= tot_test_xyz_z_only;
        }

	}
    CompositionData data; 
    data.normalized_abundances = std::move(comp_test_xyz);
    data.Xsol = tot_test_xyz_x_only;
    data.Ysol = tot_test_xyz_y_only;
    data.Zsol = tot_test_xyz_z_only;
    return data;    
    
}

double_map update_initial_composition(const CompositionData &data, double X, double Z) { 
    
    double_map comp; 
    double Hsum = 0.0;
	
	for (auto& entry : data.normalized_abundances) {


		const std::string& key = entry.first;
        double value = entry.second;

        if (key.find("H1") != std::string::npos || key.find("H2") != std::string::npos) {
            Hsum += entry.second * X; // effectively (alpha_H1 + alpha_H2)*X = alpha_H * X = X as alpha_H == 1 (this is idiot check)

        } else if (key.find("He3") != std::string::npos || key.find("He4") != std::string::npos) {
			comp[key] = entry.second * (1-X-Z); // below they did alpha_He3 * (1-X-Z) and then (1-X-Y) - comp[He3]
        } else {
			comp[key] = entry.second * Z;

        }

	}
	comp["H"] = Hsum;

	double tot=comp.sum();
	comp["Ex"] = 1 - tot;
	comp["Xsol"] = data.Xsol;
	comp["Ysol"] = data.Ysol;
	comp["Zsol"] = data.Zsol;
	global_abundance_map.Zmix = data.Zsol;

	
	global_abundance_map.comp_xa = comp; // save the entire mass fraction.
		
	return comp;

}

matrix composition_map::X() const {

	return (*this)["H"];

}

matrix composition_map::Y() const {

	return (*this)["He3"]+(*this)["He4"];

}

matrix composition_map::Z() const {

	return 1.-X()-Y();

}

double composition_map::Xsol() const {

	return (*this)["Xsol"].data()[0];
}

double composition_map::Ysol() const {

	return (*this)["Ysol"].data()[0];
}

double composition_map::Zsol() const {

	return (*this)["Zsol"].data()[0]; 
}

