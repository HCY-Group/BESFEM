#include "mfem.hpp"
#include "../include/BESFEM_All.hpp"
#include "../include/SimulationConfig.hpp"


Adjust::Adjust(Initialize_Geometry &geo, Domain_Parameters &para, const SimulationConfig &cfg)
    : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para), cfg(cfg)

{}


// Apply constant-current voltage correction for full-cell simulations.
void Adjust::AdjustConstantCurrent(double current_A, double current_C, ElectrodePotential &anode_potential, ElectrodePotential &cathode_potential,
    mfem::ParGridFunction &phA_gf, mfem::ParGridFunction &phC_gf, double &VCell)
{
    double VsrC;
    double dCrntC = std::abs(current_C - domain_parameters.gTrgI);

    // --- Current deviation scaling ---
    if (dCrntC < std::abs(domain_parameters.gTrgI) * 0.05) // if the current deviation is less than 5% of the target current
        VsrC = 0.025 * cfg.Vsr0;
    else if (dCrntC < std::abs(domain_parameters.gTrgI) * 0.10) // if the current deviation is less than 10% of the target current
        VsrC = 0.25 * cfg.Vsr0;
    else
        VsrC = 1.0 * cfg.Vsr0;

    double VsrA;
    double dCrntA = domain_parameters.gTrgI + current_A;
    
    // double dCrntA = std::abs(current_A - domain_parameters.gTrgI);
    // double dCrntA = std::abs(std::abs(current_A) - std::abs(domain_parameters.gTrgI));

    // std::cout << "dCrntA: " << dCrntA << std::endl;
    // std::cout << "5% threshold: " << std::abs(domain_parameters.gTrgI) * 0.05 << std::endl;
    // std::cout << "10% threshold: " << std::abs(domain_parameters.gTrgI) * 0.10 << std::endl;

    // --- Current deviation scaling ---
    if (dCrntA < std::abs(domain_parameters.gTrgI) * 0.05) // if the current deviation is less than 5% of the target current
        VsrA = 0.025 * cfg.Vsr0;
    else if (dCrntA < std::abs(domain_parameters.gTrgI) * 0.10) // if the current deviation is less than 10% of the target current
        VsrA = 0.25 * cfg.Vsr0;
    else
        VsrA = 1.0 * cfg.Vsr0;

    // --- Anode voltage correction ---
    double sgnA = std::copysign(1.0, domain_parameters.gTrgI - std::abs(current_A));
    double dV_A = cfg.dt * VsrA * sgnA * 0.25;
    anode_potential.AddBoundaryVoltage(dV_A);
    phA_gf += dV_A;

    // --- Cathode voltage correction ---
    double sgnC = std::copysign(1.0, domain_parameters.gTrgI - current_C);
    double dV_C = cfg.dt * VsrC * sgnC * 2.0;
    cathode_potential.AddBoundaryVoltage(-dV_C);
    phC_gf -= dV_C;

    // --- Compute overall cell voltage ---
    VCell = cathode_potential.GetBoundaryVoltage() - anode_potential.GetBoundaryVoltage();
}

