#include "Constants.hpp"


namespace Constants {

    const char* mesh_file = "Mesh_3x90_T3.mesh";
    const char* dsF_file = "dsF_3x90_T3.txt";

    // const char* mesh_file = "II_1_bin.tif";

    const double dh = 0.5e-4;                               ///< Mesh element size
    const double eps = 1.0e-6;                              ///< Small epsilon value for numerical tolerance
    const double ze = 1.0e-5;
    const double rho = 0.0501;                              ///< Lithium site density
    const double Cr = 3.0;                                  ///< C-rate for charging/discharging cycles





    const int order = 1;                                    ///< Order of the finite element basis functions



}


