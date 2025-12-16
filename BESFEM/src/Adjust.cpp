
#include "../include/Adjust.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/Domain_Parameters.hpp"
#include "mfem.hpp"
#include "../include/CnE.hpp"
#include "../include/CnA.hpp"
#include "../include/PotE.hpp"
#include "../include/PotA.hpp"
#include "../include/PotC.hpp"
#include "../include/Reaction.hpp"

Adjust::Adjust(Initialize_Geometry &geo, Domain_Parameters &para)
    : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para)

{}


// Full Cell
void Adjust::AdjustConstantCurrent(double current_A, double current_C, PotA &anode_potential, PotC &cathode_potential,
    mfem::ParGridFunction &phA_gf, mfem::ParGridFunction &phC_gf, double &VCell)
{
    double Vsr;
    double dCrnt = std::abs(current_A - domain_parameters.gTrgI);

    // --- Current deviation scaling ---
    if (dCrnt < std::abs(domain_parameters.gTrgI) * 0.05)
        Vsr = 0.025 * Constants::Vsr0;
    else if (dCrnt < std::abs(domain_parameters.gTrgI) * 0.10)
        Vsr = 0.25 * Constants::Vsr0;
    else
        Vsr = 1.0 * Constants::Vsr0;

    // --- Anode voltage correction ---
    double sgnA = std::copysign(1.0, domain_parameters.gTrgI - std::abs(current_A));
    double dV_A = Constants::dt * Vsr * sgnA * 0.10;
    anode_potential.BvA += dV_A;
    phA_gf += dV_A;

    // --- Cathode voltage correction ---
    double sgnC = std::copysign(1.0, domain_parameters.gTrgI - current_C);
    double dV_C = Constants::dt * Vsr * sgnC * 2.0;
    cathode_potential.BvC -= dV_C;
    phC_gf -= dV_C;

    // --- Compute overall cell voltage ---
    VCell = cathode_potential.BvC - anode_potential.BvA;
}

