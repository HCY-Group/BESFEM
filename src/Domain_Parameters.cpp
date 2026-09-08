#include "../include/Constants.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/Domain_Parameters.hpp"
#include "../include/MaterialProperties.hpp"
#include "../include/readtiff.h"
#include "mfem.hpp"
#include <tiffio.h>
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>
#include <sstream>

static inline void GlobalMinMax(const mfem::ParGridFunction& gf, double& gmin, double& gmax, MPI_Comm comm = MPI_COMM_WORLD)
{
    double lmin =  std::numeric_limits<double>::infinity();
    double lmax = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < gf.Size(); ++i) {
        const double v = gf(i);
        if (v < lmin) lmin = v;
        if (v > lmax) lmax = v;
    }
    MPI_Allreduce(&lmin, &gmin, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&lmax, &gmax, 1, MPI_DOUBLE, MPI_MAX, comm);
}

double gTrgI = 0.0;

Domain_Parameters::Domain_Parameters(Initialize_Geometry &geo, const SimulationConfig &cfg)
    : geometry(geo), cfg(cfg), nV(geo.nV), nE(geo.nE), nC(geo.nC), 
    pmesh(geo.parallelMesh.get()), fespace(geo.parfespace),
    particle_labels(geo.particle_labels), anode_particle_labels(geo.anode_particle_labels), cathode_particle_labels(geo.cathode_particle_labels)
{}

// Destructor
Domain_Parameters::~Domain_Parameters() {}

void Domain_Parameters::SetupDomainParameters(){

    InitializeGridFunctions();
    InterpolateDomainParameters();
    CalculatePhasePotentialsAndTargetCurrent();

    psi->SaveAsOne("psi");
    pse->SaveAsOne("pse");
    AvP->SaveAsOne("AvP");
    AvE->SaveAsOne("AvE");
    AvB->SaveAsOne("AvB");
    pmesh->SaveAsOne("pmesh");

    if (cfg.mode == sim::CellMode::FULL)
    {
        psiA->SaveAsOne("psiA");
        psiC->SaveAsOne("psiC");

        AvPA->SaveAsOne("AvPA");
        AvPC->SaveAsOne("AvPC");
    }

    PrintInfo();
}

void Domain_Parameters::InitializeGridFunctions() {

    if (!fespace) {
        throw std::runtime_error("Finite element space is not initialized.");
    }

    psi = std::make_unique<mfem::ParGridFunction>(fespace.get());
    pse = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvP = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvB = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvE = std::make_unique<mfem::ParGridFunction>(fespace.get());
    denom = std::make_unique<mfem::ParGridFunction>(fespace.get());

    if (cfg.mode == sim::CellMode::HALF)
    {
        InitializeHalfCellGridFunctions();
    }
    else 
    {
        InitializeFullCellGridFunctions();
    }    
}

void Domain_Parameters::InitializeHalfCellGridFunctions()
{
    ps.clear();
    ps.resize(particle_labels.size());

    AvPs.clear();
    AvPs.resize(particle_labels.size());

    AvEs.clear();
    AvEs.resize(particle_labels.size());

    WeightEs.clear();
    WeightEs.resize(particle_labels.size());

    for (int k = 0; k < (int)particle_labels.size(); ++k)
    {
        ps[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        AvPs[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        AvEs[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        WeightEs[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
    }

    AvP_Pairs.clear();
    AvP_Pairs.resize(particle_labels.size());
    psi_Pairs.clear();
    psi_Pairs.resize(particle_labels.size());
    WeightPairs.clear();
    WeightPairs.resize(particle_labels.size());

    for (int j = 0; j < (int)particle_labels.size(); ++j)
    {
        AvP_Pairs[j].resize(particle_labels.size());
        psi_Pairs[j].resize(particle_labels.size());
        WeightPairs[j].resize(particle_labels.size());

        for (int k = 0; k < (int)particle_labels.size(); ++k)
        {
            if (k != j)
            {
                AvP_Pairs[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
                psi_Pairs[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
                WeightPairs[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
            }
        }
    }

    tPs.clear();
    tPs.resize(particle_labels.size());

    gtPs.clear();
    gtPs.resize(particle_labels.size());

    gTrgPs.clear();
    gTrgPs.resize(particle_labels.size());

    for (int k = 0; k < (int)particle_labels.size(); ++k)
    { 
        tPs[k] = 0.0; 
        gtPs[k] = 0.0;
        gTrgPs[k] = 0.0;
    }
}

void Domain_Parameters::InitializeFullCellGridFunctions()
{
    psiA = std::make_unique<mfem::ParGridFunction>(fespace.get());
    psiC = std::make_unique<mfem::ParGridFunction>(fespace.get());

    AvPA = std::make_unique<mfem::ParGridFunction>(fespace.get());
    AvPC = std::make_unique<mfem::ParGridFunction>(fespace.get());

    denomA = std::make_unique<mfem::ParGridFunction>(fespace.get());
    denomC = std::make_unique<mfem::ParGridFunction>(fespace.get());

    const int numAnode = static_cast<int>(anode_particle_labels.size());
    const int numCathode = static_cast<int>(cathode_particle_labels.size());

    psA.resize(numAnode);
    AvPsA.resize(numAnode);
    AvEsA.resize(numAnode);
    WeightEsA.resize(numAnode);

    for (int k = 0; k < numAnode; ++k)
    {
        psA[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        AvPsA[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        AvEsA[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        WeightEsA[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
    }

    psC.resize(numCathode);
    AvPsC.resize(numCathode);
    AvEsC.resize(numCathode);
    WeightEsC.resize(numCathode);

    for (int k = 0; k < numCathode; ++k)
    {
        psC[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        AvPsC[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        AvEsC[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        WeightEsC[k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
    }

    const int num_anode_particles = static_cast<int>(psA.size());

    AvP_PairsA.clear();
    psi_PairsA.clear();
    WeightPairsA.clear();

    AvP_PairsA.resize(num_anode_particles);
    psi_PairsA.resize(num_anode_particles);

    WeightPairsA.resize(num_anode_particles);

    for (int j = 0;
        j < num_anode_particles; ++j)
    {
        AvP_PairsA[j].resize(num_anode_particles);
        psi_PairsA[j].resize(num_anode_particles);
        WeightPairsA[j].resize(num_anode_particles);

        for (int k = j + 1;
            k < num_anode_particles; ++k)
        {
            AvP_PairsA[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
            psi_PairsA[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
            WeightPairsA[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        }
    }

    const int num_cathode_particles = static_cast<int>(psC.size());

    AvP_PairsC.clear();
    psi_PairsC.clear();
    WeightPairsC.clear();

    AvP_PairsC.resize(num_cathode_particles);
    psi_PairsC.resize(num_cathode_particles);

    WeightPairsC.resize(num_cathode_particles);

    for (int j = 0;
        j < num_cathode_particles; ++j)
    {
        AvP_PairsC[j].resize(num_cathode_particles);
        psi_PairsC[j].resize(num_cathode_particles);
        WeightPairsC[j].resize(num_cathode_particles);

        for (int k = j + 1;
            k < num_cathode_particles; ++k)
        {
            AvP_PairsC[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
            psi_PairsC[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
            WeightPairsC[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
        }
    }

    tPsA.clear();
    tPsA.resize(num_anode_particles);

    tPsC.clear();
    tPsC.resize(num_cathode_particles);

    gtPsA.clear();
    gtPsA.resize(num_anode_particles);

    gtPsC.clear();
    gtPsC.resize(num_cathode_particles);

    gTrgPsA.clear();
    gTrgPsA.resize(num_anode_particles);

    gTrgPsC.clear();
    gTrgPsC.resize(num_cathode_particles);

    for (int j = 0; j < num_anode_particles; ++j)
    {
        tPsA[j] = 0.0;
        gtPsA[j] = 0.0;
        gTrgPsA[j] = 0.0;
    }

    for (int j = 0; j < num_cathode_particles; ++j)
    {
        tPsC[j] = 0.0;
        gtPsC[j] = 0.0;
        gTrgPsC[j] = 0.0;
    }

}

void Domain_Parameters::InterpolateDomainParameters() {

    nV = pmesh->GetNV();
    nE = pmesh->GetNE();
    nC = pmesh->GetElement(0)->GetNVertices();

    if (cfg.mode == sim::CellMode::HALF){
        InterpolateHalfCellMasks();
        BuildHalfCellInterfaces();
    }
    else
    {
        InterpolateFullCellMasks();
        BuildFullCellInterfaces();
    }
}

void Domain_Parameters::InterpolateHalfCellMasks()
{
    // Shared electrolyte field.
    *pse = *geometry.MaskFilterPse;
    *psi = 0.0;

    for (int k = 0; k < static_cast<int>(ps.size()); ++k)
    {
        *ps[k] = *geometry.MaskFilters[k];
        *psi += *ps[k];
    }

    // Clamp total electrode and electrolyte phase fields.
    for (int i = 0; i < psi->Size(); ++i)
    {
        (*psi)(i) = std::max(1.0e-6, std::min(1.0, (*psi)(i)));
        (*pse)(i) = std::max(1.0e-6, std::min(1.0, (*pse)(i)));
    }

    // Clamp each individual particle phase field.
    for (int k = 0; k < static_cast<int>(ps.size()); ++k)
    {
        for (int i = 0; i < ps[k]->Size(); ++i)
        {
            (*ps[k])(i) = std::max(1.0e-6, std::min(1.0, (*ps[k])(i)));
        }
    }
}

void Domain_Parameters::InterpolateFullCellMasks()
{
    // Shared electrolyte field.
    *pse = *geometry.MaskFilterPse;

    // Initialize total electrode fields.
    *psiA = 0.0;
    *psiC = 0.0;
    *psi  = 0.0;

    for (int k = 0;
         k < static_cast<int>(psA.size());
         ++k)
    {
        *psA[k] = *geometry.MaskFiltersAnode[k];
        *psiA += *psA[k];
    }

    for (int k = 0;
         k < static_cast<int>(psC.size());
         ++k)
    {
        *psC[k] = *geometry.MaskFiltersCathode[k];
        *psiC += *psC[k];
    }

    *psi = *psiA;
    *psi += *psiC;

    for (int i = 0; i < psi->Size(); ++i)
    {
        (*psiA)(i) = std::max(1.0e-6, std::min(1.0, (*psiA)(i)));
        (*psiC)(i) = std::max(1.0e-6, std::min(1.0, (*psiC)(i)));
        (*psi)(i) = std::max(1.0e-6, std::min(1.0, (*psi)(i)));
        (*pse)(i) = std::max(1.0e-6, std::min(1.0, (*pse)(i)));
    }

    for (int k = 0; k < static_cast<int>(psA.size()); ++k)
    {
        for (int i = 0; i < psA[k]->Size(); ++i)
        {
            (*psA[k])(i) = std::max(1.0e-6, std::min(1.0, (*psA[k])(i)));
        }
    }-------------------------------------------------

    for (int k = 0; k < static_cast<int>(psC.size()); ++k)
    {
        for (int i = 0; i < psC[k]->Size(); ++i)
        {
            (*psC[k])(i) = std::max(1.0e-6, std::min(1.0, (*psC[k])(i)));
        }
    }
}

void Domain_Parameters::ComputeGradientMagnitude(const mfem::ParGridFunction &phase_in, mfem::ParGridFunction &gradient_out)
{
    const int dim = pmesh->Dimension();
    gradient_out = 0.0;

    mfem::ParGridFunction derivative(fespace.get());

    for (int d = 0; d < dim; ++d)
    {
        derivative = 0.0;

        mfem::ParGridFunction phase_copy(phase_in);
        phase_copy.GetDerivative(1, d, derivative);

        for (int i = 0; i < gradient_out.Size(); ++i)
        {
            const double value = derivative(i);
            gradient_out(i) += value * value;
        }
    }

    for (int i = 0; i < gradient_out.Size(); ++i)
    {
        gradient_out(i) = std::sqrt(gradient_out(i));
    }
}

void Domain_Parameters::BuildPairInterface(mfem::ParGridFunction &out, const mfem::ParGridFunction &phase_a, const mfem::ParGridFunction &phase_b,
    const mfem::ParGridFunction &gradient_a, const mfem::ParGridFunction &gradient_b)
{
    out = phase_a;
    out *= gradient_b;

    mfem::ParGridFunction temporary(fespace.get());

    temporary = phase_b;
    temporary *= gradient_a;

    out += temporary;

    mfem::ParGridFunction overlap(fespace.get());

    overlap = phase_a;
    overlap *= phase_b;

    out *= overlap;
    out *= 4.0;

    for (int i = 0; i < out.Size(); ++i)
    {
        if (out(i) > 9000.0)
        {
            out(i) = 1.4e4;
        }
    }
}

void Domain_Parameters::BuildElectrolyteInterface(mfem::ParGridFunction &out, const mfem::ParGridFunction &electrolyte_phase, const mfem::ParGridFunction &particle_gradient)
{
    out = electrolyte_phase;
    out *= particle_gradient;
}

void Domain_Parameters::BuildPairPhaseMask(mfem::ParGridFunction &out, const mfem::ParGridFunction &phase_a, const mfem::ParGridFunction &phase_b)
{
    out = phase_a;
    out += phase_b;

    for (int i = 0; i < out.Size(); ++i)
    {
        out(i) = std::max(0.0, std::min(1.0, out(i)));
    }
}

void Domain_Parameters::ComputeInterfaceWeight(mfem::ParGridFunction &weight_out, const mfem::ParGridFunction &numerator, const mfem::ParGridFunction &denominator, const mfem::ParGridFunction *mask)
{
    weight_out = 0.0;

    const double beta = 0.8;
    const double epsilon = 1.0e-30;

    for (int i = 0; i < weight_out.Size(); ++i)
    {
        double ratio = 0.0;

        if (denominator(i) > epsilon)
        {
            ratio = numerator(i) / denominator(i);
            ratio = std::max(0.0, ratio);
        }

        weight_out(i) = std::pow(ratio, beta);
    }

    if (mask != nullptr)
    {
        weight_out *= *mask;
    }
}

void Domain_Parameters::BuildHalfCellInterfaces()
{
    ComputeGradientMagnitude(*psi, *AvP);
    ComputeGradientMagnitude(*pse, *AvE);

    for (int k = 0; k < static_cast<int>(ps.size()); ++k)
    {
        ComputeGradientMagnitude(*ps[k], *AvPs[k]);
    }

    for (int j = 0; j < static_cast<int>(ps.size()); ++j)
    {
        for (int k = j + 1; k < static_cast<int>(ps.size()); ++k)
        {
            BuildPairInterface(*AvP_Pairs[j][k], *ps[j], *ps[k], *AvPs[j], *AvPs[k]);
            BuildPairPhaseMask(*psi_Pairs[j][k], *ps[j], *ps[k]);

            std::ostringstream filename;
            filename << "AvP_Pair_" << j << "_" << k;
            AvP_Pairs[j][k]->SaveAsOne(filename.str().c_str());
        }
    }

    for (int k = 0; k < static_cast<int>(ps.size()); ++k)
    {
        BuildElectrolyteInterface(*AvEs[k], *pse, *AvPs[k]);
    }

    *denom = 0.0;

    for (int j = 0; j < static_cast<int>(ps.size()); ++j)
    {
        for (int k = j + 1; k < static_cast<int>(ps.size()); ++k)
        {
            *denom += *AvP_Pairs[j][k];
        }
    }

    for (int k = 0; k < static_cast<int>(ps.size()); ++k)
    {
        *denom += *AvEs[k];
    }

    for (int k = 0; k < static_cast<int>(ps.size()); ++k)
    {
        ComputeInterfaceWeight(*WeightEs[k], *AvEs[k], *denom);
    }

    for (int j = 0; j < static_cast<int>(ps.size()); ++j)
    {
        for (int k = j + 1; k < static_cast<int>(ps.size()); ++k)
        {
            ComputeInterfaceWeight(*WeightPairs[j][k], *AvP_Pairs[j][k], *denom, psi_Pairs[j][k].get());
        }
    }
}

void Domain_Parameters::BuildFullCellInterfaces()
{
    ComputeGradientMagnitude(*psi, *AvP);
    ComputeGradientMagnitude(*psiA, *AvPA);
    ComputeGradientMagnitude(*psiC, *AvPC);
    ComputeGradientMagnitude(*pse, *AvE);

    for (int k = 0; k < static_cast<int>(psA.size()); ++k)
    {
        ComputeGradientMagnitude(*psA[k], *AvPsA[k]);
    }

    for (int k = 0; k < static_cast<int>(psC.size()); ++k)
    {
        ComputeGradientMagnitude(*psC[k], *AvPsC[k]);
    }

    for (int j = 0; j < static_cast<int>(psA.size()); ++j)
    {
        for (int k = j + 1; k < static_cast<int>(psA.size()); ++k)
        {
            BuildPairInterface(*AvP_PairsA[j][k], *psA[j], *psA[k], *AvPsA[j], *AvPsA[k]);
            BuildPairPhaseMask(*psi_PairsA[j][k], *psA[j], *psA[k]);
        }
    }

    for (int j = 0; j < static_cast<int>(psC.size()); ++j)
    {
        for (int k = j + 1; k < static_cast<int>(psC.size()); ++k)
        {
            BuildPairInterface(*AvP_PairsC[j][k], *psC[j], *psC[k], *AvPsC[j], *AvPsC[k]);
            BuildPairPhaseMask(*psi_PairsC[j][k], *psC[j], *psC[k]);
        }
    }

    for (int k = 0; k < static_cast<int>(psA.size()); ++k)
    {
        BuildElectrolyteInterface(*AvEsA[k], *pse, *AvPsA[k]);
    }

    for (int k = 0; k < static_cast<int>(psC.size()); ++k)
    {
        BuildElectrolyteInterface(*AvEsC[k], *pse, *AvPsC[k]);
    }

    *denomA = 0.0;

    for (int j = 0; j < static_cast<int>(psA.size()); ++j)
    {
        for (int k = j + 1; k < static_cast<int>(psA.size()); ++k)
        {
            *denomA += *AvP_PairsA[j][k];
        }
    }

    for (int k = 0; k < static_cast<int>(psA.size()); ++k)
    {
        *denomA += *AvEsA[k];
    }

    *denomC = 0.0;

    for (int j = 0; j < static_cast<int>(psC.size()); ++j)
    {
        for (int k = j + 1; k < static_cast<int>(psC.size()); ++k)
        {
            *denomC += *AvP_PairsC[j][k];
        }
    }

    for (int k = 0; k < static_cast<int>(psC.size()); ++k)
    {
        *denomC += *AvEsC[k];
    }

    for (int k = 0; k < static_cast<int>(psA.size()); ++k)
    {
        ComputeInterfaceWeight(*WeightEsA[k], *AvEsA[k], *denomA);
    }

    for (int k = 0; k < static_cast<int>(psC.size()); ++k)
    {
        ComputeInterfaceWeight(*WeightEsC[k], *AvEsC[k], *denomC);
    }

    for (int j = 0; j < static_cast<int>(psA.size()); ++j)
    {
        for (int k = j + 1; k < static_cast<int>(psA.size()); ++k)
        {
            ComputeInterfaceWeight(*WeightPairsA[j][k], *AvP_PairsA[j][k], *denomA, psi_PairsA[j][k].get());
        }
    }

    for (int j = 0; j < static_cast<int>(psC.size()); ++j)
    {
        for (int k = j + 1; k < static_cast<int>(psC.size()); ++k)
        {
            ComputeInterfaceWeight(*WeightPairsC[j][k], *AvP_PairsC[j][k], *denomC, psi_PairsC[j][k].get());
        }
    }
}


void Domain_Parameters::CalculateTotals(const mfem::ParGridFunction &grid_function, const mfem::Vector &element_volumes, double &local_total, double &global_total)
{
    local_total = 0.0;

    for (int ei = 0; ei < pmesh->GetNE(); ++ei)
    {
        mfem::Array<double> nodal_values;
        grid_function.GetNodalValues(ei, nodal_values);

        if (nodal_values.Size() == 0)
        {
            continue;
        }

        double average_value = 0.0;

        for (int j = 0; j < nodal_values.Size(); ++j)
        {
            average_value += nodal_values[j];
        }

        average_value /= static_cast<double>(nodal_values.Size());
        local_total += average_value * element_volumes(ei);
    }

    MPI_Allreduce(&local_total, &global_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void Domain_Parameters::CalculateTotalPhaseField(const mfem::ParGridFunction &grid_function, double &local_total, double &global_total)
{
    const int local_element_count = pmesh->GetNE();
    EVol.SetSize(local_element_count);

    for (int ei = 0; ei < local_element_count; ++ei)
    {
        EVol(ei) = pmesh->GetElementVolume(ei);
    }

    CalculateTotals(grid_function, EVol, local_total, global_total);
}

void Domain_Parameters::CalculatePhasePotentialsAndTargetCurrent()
{
    // Electrolyte is shared by both half-cell and full-cell modes.
    CalculateTotalPhaseField(*pse, tPse, gtPse);

    if (cfg.mode == sim::CellMode::HALF)
    {
        CalculateHalfCellPhasePotentialsAndTargetCurrent();
    }
    else
    {
        CalculateFullCellPhasePotentialsAndTargetCurrent();
    }
}

void Domain_Parameters::CalculateHalfCellPhasePotentialsAndTargetCurrent()
{
    const std::vector<sim::MaterialType> &active_materials = (cfg.half_electrode == sim::Electrode::CATHODE) ? cfg.cathode_materials : cfg.anode_materials;

    MFEM_VERIFY(active_materials.size() == ps.size(), "Half-cell material count does not match particle count.");

    gTrgI = 0.0;
    CalculateTotalPhaseField(*psi, tPsi, gtPsi);

    for (std::size_t k = 0; k < ps.size(); ++k)
    {
        CalculateTotalPhaseField(*ps[k], tPs[k], gtPs[k]);
        CalculateTargetCurrent(tPs[k], gTrgPs[k], active_materials[k]);
        gTrgI += gTrgPs[k];
    }
}

void Domain_Parameters::CalculateFullCellPhasePotentialsAndTargetCurrent()
{
    gTrgIA = 0.0;
    gTrgIC = 0.0;
    gTrgI = 0.0;

    CalculateTotalPhaseField(*psiA, tPsiA,gtPsiA);
    CalculateTotalPhaseField(*psiC, tPsiC, gtPsiC);

    CalculateTotalPhaseField(*psi, tPsi, gtPsi);

    for (std::size_t k = 0; k < psA.size(); ++k)
    {
        CalculateTotalPhaseField(*psA[k], tPsA[k], gtPsA[k]);
        CalculateTargetCurrent(tPsA[k], gTrgPsA[k], cfg.anode_materials[k]);
        gTrgIA += gTrgPsA[k];
    }

    for (std::size_t k = 0; k < psC.size(); ++k)
    {
        CalculateTotalPhaseField(*psC[k], tPsC[k], gtPsC[k]);
        CalculateTargetCurrent(tPsC[k], gTrgPsC[k], cfg.cathode_materials[k]);
        gTrgIC += gTrgPsC[k];
    }

    gTrgI = gTrgIC; // use cathode target current in the full cell
}

void Domain_Parameters::CalculateTargetCurrent(double local_phase_volume, double &global_target_current, sim::MaterialType material)
{
    const double rho = MaterialProperties::SiteDensity(material);
    const double local_target_current = local_phase_volume * rho * (0.95 - 0.3) / (3600.0 / cfg.Cr);
    MPI_Allreduce(&local_target_current, &global_target_current, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void Domain_Parameters::PrintInfo()
{
    if (mfem::Mpi::WorldRank() != 0)
    {
        return;
    }

    std::cout << "Total solid phase: " << gtPsi << '\n' << "Total electrolyte phase: " << gtPse << '\n';

    if (cfg.mode == sim::CellMode::HALF)
    {
        std::cout << "Target Current: " << gTrgI << '\n';

        for (std::size_t k = 0; k < gtPs.size(); ++k)
        {
            std::cout << "Particle " << k << " phase total: " << gtPs[k] << ", target current: " << gTrgPs[k] << '\n';
        }
    }
    else
    {
        std::cout << "Total anode phase: " << gtPsiA << '\n' << "Total cathode phase: " << gtPsiC << '\n'
            << "Anode capacity-based current: " << gTrgIA << '\n' << "Cathode capacity-based current: " << gTrgIC << '\n'
            << "Selected full-cell target current: " << gTrgI << '\n';

        for (std::size_t k = 0; k < gtPsA.size(); ++k)
        {
            std::cout << "Anode particle " << k << " phase total: "
                << gtPsA[k] << ", target current: " << gTrgPsA[k] << '\n';
        }

        for (std::size_t k = 0; k < gtPsC.size(); ++k)
        {
            std::cout << "Cathode particle " << k  << " phase total: " << gtPsC[k] << ", target current: " << gTrgPsC[k] << '\n';
        }
    }
}

