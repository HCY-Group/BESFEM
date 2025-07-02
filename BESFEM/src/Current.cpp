/**
 * @file Current.cpp
 * @brief Implementation of the Current class for controlling current in battery simulations.
 */

#include "../include/Current.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/Domain_Parameters.hpp"
#include "mfem.hpp"
#include "../include/CnE.hpp"
#include "../include/CnP.hpp"
#include "../include/PotE.hpp"
#include "../include/PotP.hpp"
#include "../include/Reaction.hpp"

Current::Current(Initialize_Geometry &geo, Domain_Parameters &para)
    : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para)

{}

void Current::Constant(mfem::ParGridFunction &phx, double &global_current){

    sgn = copysign(1, gTrgI - global_current);     // Compute the sign of the difference between the target and global current
    dV = Constants::dt * Constants::Vsr * sgn;     // Compute the voltage adjustment step
    
    // Adjust the boundary value and particle potential
    BvP -= dV; 
    phx -= dV;

    
}


