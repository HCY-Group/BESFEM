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

static inline void GlobalMinMax(const mfem::ParGridFunction& gf,
                                double& gmin, double& gmax,
                                MPI_Comm comm = MPI_COMM_WORLD)
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

        // denomA->SaveAsOne("denomA");
        // denomC->SaveAsOne("denomC");

        // for (int k = 0;
        //     k < static_cast<int>(AvEsA.size());
        //     ++k)
        // {
        //     std::ostringstream filename;
        //     filename << "AvEsA_" << k;
        //     AvEsA[k]->SaveAsOne(filename.str().c_str());
        // }

        // for (int k = 0;
        //     k < static_cast<int>(AvEsC.size());
        //     ++k)
        // {
        //     std::ostringstream filename;
        //     filename << "AvEsC_" << k;
        //     AvEsC[k]->SaveAsOne(filename.str().c_str());
        // }
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
        j < num_anode_particles;
        ++j)
    {
        AvP_PairsA[j].resize(num_anode_particles);
        psi_PairsA[j].resize(num_anode_particles);
        WeightPairsA[j].resize(num_anode_particles);

        for (int k = j + 1;
            k < num_anode_particles;
            ++k)
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
        j < num_cathode_particles;
        ++j)
    {
        AvP_PairsC[j].resize(num_cathode_particles);
        psi_PairsC[j].resize(num_cathode_particles);
        WeightPairsC[j].resize(num_cathode_particles);

        for (int k = j + 1;
            k < num_cathode_particles;
            ++k)
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

    // nV = pmesh->GetNV();

    MFEM_VERIFY(pmesh, "Parallel mesh is not initialized.");
    MFEM_VERIFY(geometry.MaskFilterPse, "Electrolyte mask is not initialized.");
    MFEM_VERIFY(static_cast<int>(geometry.MaskFilters.size()) == static_cast<int>(ps.size()),
                "Particle mask count does not match particle field count.");

    nV = pmesh->GetNV();
    nE = pmesh->GetNE();
    nC = pmesh->GetElement(0)->GetNVertices();

    if (cfg.mode == sim::CellMode::HALF){
        MFEM_VERIFY(geometry.MaskFilters.size() == ps.size(), "Half-cell particle mask count does not match the particle field count.");
        InterpolateHalfCellMasks();
        BuildHalfCellInterfaces();

    }
    else
    {
        MFEM_VERIFY(geometry.MaskFiltersAnode.size() == psA.size(), "Full-cell anode particle mask count does not match the anode particle field count.");
        MFEM_VERIFY(geometry.MaskFiltersCathode.size() == psC.size(), "Full-cell cathode particle mask count does not match the cathode particle field count.");
        InterpolateFullCellMasks();
        BuildFullCellInterfaces();
    }

    // ApplyAMR();

    // if (cfg.mode == sim::CellMode::HALF){
    //     BuildHalfCellInterfaces();
    // }
    // else
    // {
    //     BuildFullCellInterfaces();
    // }
}

void Domain_Parameters::InterpolateHalfCellMasks()
{
    MFEM_VERIFY(geometry.MaskFilterPse, "Half-cell electrolyte mask is not initialized.");

    MFEM_VERIFY(static_cast<int>(geometry.MaskFilters.size()) == static_cast<int>(ps.size()),
        "Half-cell particle mask count does not match the particle field count.");

    // Shared electrolyte field.
    *pse = *geometry.MaskFilterPse;
    *psi = 0.0;

    for (int k = 0;
         k < static_cast<int>(ps.size());
         ++k)
    {
        MFEM_VERIFY(
            geometry.MaskFilters[k],
            "A half-cell particle mask is not initialized.");

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
    for (int k = 0;
         k < static_cast<int>(ps.size());
         ++k)
    {
        for (int i = 0;
             i < ps[k]->Size();
             ++i)
        {
            (*ps[k])(i) = std::max(1.0e-6, std::min(1.0, (*ps[k])(i)));
        }
    }
}

void Domain_Parameters::InterpolateFullCellMasks()
{
    MFEM_VERIFY(geometry.MaskFilterPse, "Full-cell electrolyte mask is not initialized.");
    MFEM_VERIFY(geometry.MaskFilterAnode, "Full-cell anode mask is not initialized.");
    MFEM_VERIFY(geometry.MaskFilterCathode, "Full-cell cathode mask is not initialized.");

    MFEM_VERIFY(static_cast<int>(geometry.MaskFiltersAnode.size()) == static_cast<int>(psA.size()),
        "Anode particle-mask count does not match the anode particle-field count.");

    MFEM_VERIFY(static_cast<int>(geometry.MaskFiltersCathode.size()) == static_cast<int>(psC.size()),
        "Cathode particle-mask count does not match the cathode particle-field count.");

    // Shared electrolyte field.
    *pse = *geometry.MaskFilterPse;

    // Initialize total electrode fields.
    *psiA = 0.0;
    *psiC = 0.0;
    *psi  = 0.0;

    // -------------------------------------------------
    // Anode particle fields
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(psA.size());
         ++k)
    {
        MFEM_VERIFY(geometry.MaskFiltersAnode[k], "An anode particle mask is not initialized.");
        *psA[k] = *geometry.MaskFiltersAnode[k];
        *psiA += *psA[k];
    }

    // -------------------------------------------------
    // Cathode particle fields
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(psC.size());
         ++k)
    {
        MFEM_VERIFY(geometry.MaskFiltersCathode[k], "A cathode particle mask is not initialized.");
        *psC[k] = *geometry.MaskFiltersCathode[k];
        *psiC += *psC[k];
    }

    // Total solid field used for shared calculations
    // such as AMR and total solid volume.
    *psi = *psiA;
    *psi += *psiC;

    // -------------------------------------------------
    // Clamp total phase fields
    // -------------------------------------------------

    for (int i = 0;
         i < psi->Size();
         ++i)
    {
        (*psiA)(i) = std::max(1.0e-6, std::min(1.0, (*psiA)(i)));
        (*psiC)(i) = std::max(1.0e-6, std::min(1.0, (*psiC)(i)));
        (*psi)(i) = std::max(1.0e-6, std::min(1.0, (*psi)(i)));
        (*pse)(i) = std::max(1.0e-6, std::min(1.0, (*pse)(i)));
    }

    // -------------------------------------------------
    // Clamp individual anode particle fields
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(psA.size());
         ++k)
    {
        for (int i = 0;
             i < psA[k]->Size();
             ++i)
        {
            (*psA[k])(i) = std::max(1.0e-6, std::min(1.0, (*psA[k])(i)));
        }
    }

    // -------------------------------------------------
    // Clamp individual cathode particle fields
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(psC.size());
         ++k)
    {
        for (int i = 0;
             i < psC[k]->Size();
             ++i)
        {
            (*psC[k])(i) = std::max(1.0e-6, std::min(1.0, (*psC[k])(i)));
        }
    }
}

// void Domain_Parameters::ApplyAMR()
// {
//     if (cfg.amr_levels <= 0)
//     {
//         return;
//     }

//     MFEM_VERIFY(pmesh, "Parallel mesh is not initialized.");
//     MFEM_VERIFY(fespace, "Finite element space is not initialized.");
//     MFEM_VERIFY(psi, "Total solid phase field is not initialized.");

//     const double outer_half_width = 0.45;

//     for (int lev = 0;
//          lev < cfg.amr_levels;
//          ++lev)
//     {
//         mfem::Array<int> refinement_list;

//         const double band_fraction = static_cast<double>(cfg.amr_levels - lev) / static_cast<double>(cfg.amr_levels);
//         const double half_width = outer_half_width * band_fraction;
//         const double psi_lower = 0.5 - half_width;
//         const double psi_upper = 0.5 + half_width;

//         // -------------------------------------------------
//         // Mark elements using the total solid phase field
//         // -------------------------------------------------

//         for (int ei = 0;
//              ei < pmesh->GetNE();
//              ++ei)
//         {
//             mfem::Array<double> psi_values;

//             psi->GetNodalValues(ei, psi_values);
//             double psi_average = 0.0;

//             for (int j = 0;
//                  j < psi_values.Size();
//                  ++j)
//             {
//                 psi_average += psi_values[j];
//             }

//             if (psi_values.Size() > 0)
//             {
//                 psi_average /= static_cast<double>(psi_values.Size());
//             }

//             if (psi_average > psi_lower && psi_average < psi_upper)
//             {
//                 refinement_list.Append(ei);
//             }
//         }

//         // -------------------------------------------------
//         // Global element counts
//         // -------------------------------------------------

//         const int local_marked = refinement_list.Size();

//         const int local_elements = pmesh->GetNE();

//         int global_marked = 0;
//         int global_elements = 0;

//         MPI_Allreduce(&local_marked, &global_marked, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//         MPI_Allreduce(&local_elements, &global_elements, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

//         if (mfem::Mpi::WorldRank() == 0)
//         {
//             std::cout << "[AMR] band " << lev + 1
//                 << ": psi range = (" << psi_lower << ", " << psi_upper << "), marked "
//                 << global_marked << " / " << global_elements << " elements globally"
//                 << std::endl;
//         }

//         if (global_marked == 0)
//         {
//             if (mfem::Mpi::WorldRank() == 0)
//             {
//                 std::cout << "[AMR] No elements found in band " << lev + 1
//                     << ". Stopping refinement." << std::endl;
//             }

//             break;
//         }

//         // -------------------------------------------------
//         // Refine mesh and update the finite element space
//         // -------------------------------------------------

//         pmesh->GeneralRefinement(refinement_list, 1);
//         fespace->Update();

//         // -------------------------------------------------
//         // Update grid functions
//         // -------------------------------------------------

//         auto UpdateGridFunction = [](std::unique_ptr<mfem::ParGridFunction> &field)
//         {
//             if (field)
//             {
//                 field->Update();
//             }
//         };

//         auto UpdateGridFunctionVector = [&](std::vector<std::unique_ptr<mfem::ParGridFunction>> &fields)
//         {
//             for (auto &field : fields)
//             {
//                 UpdateGridFunction(field);
//             }
//         };

//         auto UpdateGridFunctionMatrix =[&](std::vector<std::vector<std::unique_ptr<mfem::ParGridFunction>>> &fields)
//         {
//             for (auto &row : fields)
//             {
//                 for (auto &field : row)
//                 {
//                     UpdateGridFunction(field);
//                 }
//             }
//         };

//         // Shared fields.
//         UpdateGridFunction(psi);
//         UpdateGridFunction(pse);

//         UpdateGridFunction(AvP);
//         UpdateGridFunction(AvB);
//         UpdateGridFunction(AvE);
//         UpdateGridFunction(denom);

//         // -------------------------------------------------
//         // Half-cell fields
//         // -------------------------------------------------

//         if (cfg.mode == sim::CellMode::HALF)
//         {
//             UpdateGridFunctionVector(ps);
//             UpdateGridFunctionVector(AvPs);
//             UpdateGridFunctionVector(AvEs);
//             UpdateGridFunctionVector(WeightEs);

//             UpdateGridFunctionMatrix(AvP_Pairs);
//             UpdateGridFunctionMatrix(psi_Pairs);
//             UpdateGridFunctionMatrix(WeightPairs);
//         }

//         // -------------------------------------------------
//         // Full-cell fields
//         // -------------------------------------------------

//         else
//         {
//             UpdateGridFunction(psiA);
//             UpdateGridFunction(psiC);

//             UpdateGridFunctionVector(psA);
//             UpdateGridFunctionVector(psC);

//             UpdateGridFunctionVector(AvPsA);
//             UpdateGridFunctionVector(AvPsC);

//             UpdateGridFunctionVector(AvEsA);
//             UpdateGridFunctionVector(AvEsC);

//             UpdateGridFunctionVector(WeightEsA);
//             UpdateGridFunctionVector(WeightEsC);

//             UpdateGridFunctionMatrix(AvP_PairsA);
//             UpdateGridFunctionMatrix(AvP_PairsC);

//             UpdateGridFunctionMatrix(psi_PairsA);
//             UpdateGridFunctionMatrix(psi_PairsC);

//             UpdateGridFunctionMatrix(WeightPairsA);
//             UpdateGridFunctionMatrix(WeightPairsC);

//             UpdateGridFunction(AvPA);
//             UpdateGridFunction(AvPC);

//             UpdateGridFunction(denomA);
//             UpdateGridFunction(denomC);
//         }

//         // -------------------------------------------------
//         // Update geometry masks attached to the same space
//         // -------------------------------------------------

//         UpdateGridFunction(geometry.MaskFilterPse);

//         if (cfg.mode == sim::CellMode::HALF)
//         {
//             UpdateGridFunction(geometry.MaskFilter);

//             for (auto &field : geometry.MaskFilters)
//             {
//                 UpdateGridFunction(field);
//             }
//         }
//         else
//         {
//             UpdateGridFunction(geometry.MaskFilterAnode);

//             UpdateGridFunction(geometry.MaskFilterCathode);

//             for (auto &field : geometry.MaskFiltersAnode)
//             {
//                 UpdateGridFunction(field);
//             }

//             for (auto &field : geometry.MaskFiltersCathode)
//             {
//                 UpdateGridFunction(field);
//             }
//         }

//         // -------------------------------------------------
//         // Update stored mesh dimensions
//         // -------------------------------------------------

//         nV = pmesh->GetNV();
//         nE = pmesh->GetNE();

//         if (nE > 0)
//         {
//             nC = pmesh->GetElement(0)->GetNVertices();
//         }

//         geometry.nV = nV;
//         geometry.nE = nE;
//         geometry.nC = nC;

//         // -------------------------------------------------
//         // Determine element-size range
//         // -------------------------------------------------

//         double local_hmin = std::numeric_limits<double>::max();

//         double local_hmax = 0.0;

//         for (int ei = 0;
//              ei < pmesh->GetNE();
//              ++ei)
//         {
//             const double element_size = pmesh->GetElementSize(ei);

//             local_hmin = std::min(local_hmin, element_size);

//             local_hmax = std::max(local_hmax, element_size);
//         }

//         double global_hmin = 0.0;
//         double global_hmax = 0.0;

//         MPI_Allreduce(&local_hmin, &global_hmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//         MPI_Allreduce(&local_hmax, &global_hmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

//         const int local_elements_after = pmesh->GetNE();
//         int global_elements_after = 0;

//         MPI_Allreduce(&local_elements_after, &global_elements_after, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

//         if (mfem::Mpi::WorldRank() == 0)
//         {
//             std::cout << "[AMR] band " << lev + 1 << " element-size range after refinement: "
//                 << global_hmin << " to " << global_hmax << std::endl;

//             std::cout << "[AMR] band " << lev + 1 << " complete: " << global_elements_after << " total elements" << std::endl;
//         }
//     }
// }

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

        for (int i = 0;
             i < gradient_out.Size();
             ++i)
        {
            const double value = derivative(i);
            gradient_out(i) += value * value;
        }
    }

    for (int i = 0;
         i < gradient_out.Size();
         ++i)
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

    for (int i = 0;
         i < out.Size();
         ++i)
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

    for (int i = 0;
         i < out.Size();
         ++i)
    {
        out(i) = std::max(0.0, std::min(1.0, out(i)));
    }
}

void Domain_Parameters::ComputeInterfaceWeight(mfem::ParGridFunction &weight_out, const mfem::ParGridFunction &numerator, const mfem::ParGridFunction &denominator, const mfem::ParGridFunction *mask)
{
    weight_out = 0.0;

    const double beta = 0.8;
    const double epsilon = 1.0e-30;

    MFEM_VERIFY(weight_out.Size() == numerator.Size(), "Weight and numerator sizes do not match.");
    MFEM_VERIFY(denominator.Size() == numerator.Size(), "Denominator and numerator sizes do not match.");

    for (int i = 0;
         i < weight_out.Size();
         ++i)
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
    MFEM_VERIFY(psi, "Half-cell total solid field is not initialized.");
    MFEM_VERIFY(pse, "Half-cell electrolyte field is not initialized.");
    MFEM_VERIFY(AvP, "Half-cell total solid interface field is not initialized.");
    MFEM_VERIFY(AvE, "Half-cell electrolyte interface field is not initialized.");
    MFEM_VERIFY(denom, "Half-cell interface denominator is not initialized.");
    MFEM_VERIFY(ps.size() == AvPs.size(), "Half-cell particle and gradient-field counts do not match.");
    MFEM_VERIFY(ps.size() == AvEs.size(), "Half-cell particle and electrolyte-interface counts do not match.");
    MFEM_VERIFY(ps.size() == WeightEs.size(), "Half-cell particle and electrolyte-weight counts do not match.");

    // -------------------------------------------------
    // Total solid and electrolyte gradient magnitudes
    // -------------------------------------------------

    ComputeGradientMagnitude(*psi, *AvP);
    ComputeGradientMagnitude(*pse, *AvE);

    // -------------------------------------------------
    // Individual particle gradient magnitudes
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(ps.size());
         ++k)
    {
        ComputeGradientMagnitude(*ps[k], *AvPs[k]);
    }

    // -------------------------------------------------
    // Particle-particle interfaces
    // -------------------------------------------------

    for (int j = 0;
         j < static_cast<int>(ps.size());
         ++j)
    {
        for (int k = j + 1;
             k < static_cast<int>(ps.size());
             ++k)
        {
            MFEM_VERIFY(AvP_Pairs[j][k], "Half-cell particle-pair interface field is missing.");
            MFEM_VERIFY(psi_Pairs[j][k], "Half-cell particle-pair phase field is missing.");

            BuildPairInterface(*AvP_Pairs[j][k], *ps[j], *ps[k], *AvPs[j], *AvPs[k]);
            BuildPairPhaseMask(*psi_Pairs[j][k], *ps[j], *ps[k]);

            std::ostringstream filename;
            filename << "AvP_Pair_" << j << "_" << k;
            AvP_Pairs[j][k]->SaveAsOne(filename.str().c_str());
            
        }
    }

    // -------------------------------------------------
    // Particle-electrolyte interfaces
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(ps.size());
         ++k)
    {
        BuildElectrolyteInterface(*AvEs[k], *pse, *AvPs[k]);
    }

    // -------------------------------------------------
    // Build common half-cell denominator
    // -------------------------------------------------

    *denom = 0.0;

    for (int j = 0;
         j < static_cast<int>(ps.size());
         ++j)
    {
        for (int k = j + 1;
             k < static_cast<int>(ps.size());
             ++k)
        {
            *denom += *AvP_Pairs[j][k];
        }
    }

    for (int k = 0;
         k < static_cast<int>(ps.size());
         ++k)
    {
        *denom += *AvEs[k];
    }

    // -------------------------------------------------
    // Particle-electrolyte weights
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(ps.size());
         ++k)
    {
        ComputeInterfaceWeight(*WeightEs[k], *AvEs[k], *denom);
    }

    // -------------------------------------------------
    // Particle-particle weights
    // -------------------------------------------------

    for (int j = 0;
         j < static_cast<int>(ps.size());
         ++j)
    {
        for (int k = j + 1;
             k < static_cast<int>(ps.size());
             ++k)
        {
            MFEM_VERIFY(WeightPairs[j][k], "Half-cell particle-pair weight field is missing.");
            ComputeInterfaceWeight(*WeightPairs[j][k], *AvP_Pairs[j][k], *denom, psi_Pairs[j][k].get());
            // std::ostringstream filename;
            // filename << "WeightPair_" << j << "_" << k;
            // WeightPairs[j][k]->SaveAsOne(filename.str().c_str());
        }
    }
}

void Domain_Parameters::BuildFullCellInterfaces()
{
    MFEM_VERIFY(psi, "Full-cell total solid field is not initialized.");
    MFEM_VERIFY(psiA, "Full-cell anode phase field is not initialized.");
    MFEM_VERIFY(psiC, "Full-cell cathode phase field is not initialized.");
    MFEM_VERIFY(pse, "Full-cell electrolyte phase field is not initialized.");
    MFEM_VERIFY(AvP, "Full-cell total solid interface field is not initialized.");
    MFEM_VERIFY(AvPA, "Full-cell anode interface field is not initialized.");
    MFEM_VERIFY(AvPC, "Full-cell cathode interface field is not initialized.");
    MFEM_VERIFY(AvE, "Full-cell electrolyte interface field is not initialized.");
    MFEM_VERIFY(denomA, "Full-cell anode interface denominator is not initialized.");
    MFEM_VERIFY(denomC, "Full-cell cathode interface denominator is not initialized.");
    MFEM_VERIFY(psA.size() == AvPsA.size(), "Anode particle and gradient-field counts do not match.");
    MFEM_VERIFY(psA.size() == AvEsA.size(), "Anode particle and electrolyte-interface counts do not match.");
    MFEM_VERIFY(psA.size() == WeightEsA.size(), "Anode particle and electrolyte-weight counts do not match.");
    MFEM_VERIFY(psC.size() == AvPsC.size(), "Cathode particle and gradient-field counts do not match.");
    MFEM_VERIFY(psC.size() == AvEsC.size(), "Cathode particle and electrolyte-interface counts do not match.");
    MFEM_VERIFY(psC.size() == WeightEsC.size(), "Cathode particle and electrolyte-weight counts do not match.");

    // -------------------------------------------------
    // Total phase-field gradient magnitudes
    // -------------------------------------------------

    ComputeGradientMagnitude(*psi, *AvP);
    ComputeGradientMagnitude(*psiA, *AvPA);
    ComputeGradientMagnitude(*psiC, *AvPC);
    ComputeGradientMagnitude(*pse, *AvE);

    // -------------------------------------------------
    // Anode particle gradient magnitudes
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(psA.size());
         ++k)
    {
        ComputeGradientMagnitude(*psA[k], *AvPsA[k]);
    }

    // -------------------------------------------------
    // Cathode particle gradient magnitudes
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(psC.size());
         ++k)
    {
        ComputeGradientMagnitude(*psC[k], *AvPsC[k]);
    }

    // -------------------------------------------------
    // Anode particle-particle interfaces
    // -------------------------------------------------

    for (int j = 0;
         j < static_cast<int>(psA.size());
         ++j)
    {
        for (int k = j + 1;
             k < static_cast<int>(psA.size());
             ++k)
        {
            MFEM_VERIFY(AvP_PairsA[j][k], "Anode particle-pair interface field is missing.");
            MFEM_VERIFY(psi_PairsA[j][k], "Anode particle-pair phase field is missing.");

            BuildPairInterface(*AvP_PairsA[j][k], *psA[j], *psA[k], *AvPsA[j], *AvPsA[k]);
            BuildPairPhaseMask(*psi_PairsA[j][k], *psA[j], *psA[k]);
        }
    }

    // -------------------------------------------------
    // Cathode particle-particle interfaces
    // -------------------------------------------------

    for (int j = 0;
         j < static_cast<int>(psC.size());
         ++j)
    {
        for (int k = j + 1;
             k < static_cast<int>(psC.size());
             ++k)
        {
            MFEM_VERIFY(AvP_PairsC[j][k], "Cathode particle-pair interface field is missing.");
            MFEM_VERIFY(psi_PairsC[j][k], "Cathode particle-pair phase field is missing.");

            BuildPairInterface(*AvP_PairsC[j][k], *psC[j], *psC[k], *AvPsC[j], *AvPsC[k]);
            BuildPairPhaseMask(*psi_PairsC[j][k], *psC[j], *psC[k]);
        }
    }

    // -------------------------------------------------
    // Anode-electrolyte interfaces
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(psA.size());
         ++k)
    {
        BuildElectrolyteInterface(*AvEsA[k], *pse, *AvPsA[k]);
    }

    // -------------------------------------------------
    // Cathode-electrolyte interfaces
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(psC.size());
         ++k)
    {
        BuildElectrolyteInterface(*AvEsC[k], *pse, *AvPsC[k]);
    }

    // -------------------------------------------------
    // Anode interface denominator
    // -------------------------------------------------

    *denomA = 0.0;

    for (int j = 0;
         j < static_cast<int>(psA.size());
         ++j)
    {
        for (int k = j + 1;
             k < static_cast<int>(psA.size());
             ++k)
        {
            *denomA += *AvP_PairsA[j][k];
        }
    }

    for (int k = 0;
         k < static_cast<int>(psA.size());
         ++k)
    {
        *denomA += *AvEsA[k];
    }

    // -------------------------------------------------
    // Cathode interface denominator
    // -------------------------------------------------

    *denomC = 0.0;

    for (int j = 0;
         j < static_cast<int>(psC.size());
         ++j)
    {
        for (int k = j + 1;
             k < static_cast<int>(psC.size());
             ++k)
        {
            *denomC += *AvP_PairsC[j][k];
        }
    }

    for (int k = 0;
         k < static_cast<int>(psC.size());
         ++k)
    {
        *denomC += *AvEsC[k];
    }

    // -------------------------------------------------
    // Anode-electrolyte weights
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(psA.size());
         ++k)
    {
        ComputeInterfaceWeight(*WeightEsA[k], *AvEsA[k], *denomA);
    }

    // -------------------------------------------------
    // Cathode-electrolyte weights
    // -------------------------------------------------

    for (int k = 0;
         k < static_cast<int>(psC.size());
         ++k)
    {
        ComputeInterfaceWeight(*WeightEsC[k], *AvEsC[k], *denomC);
    }

    // -------------------------------------------------
    // Anode particle-particle weights
    // -------------------------------------------------

    for (int j = 0;
         j < static_cast<int>(psA.size());
         ++j)
    {
        for (int k = j + 1;
             k < static_cast<int>(psA.size());
             ++k)
        {
            MFEM_VERIFY(WeightPairsA[j][k], "Anode particle-pair weight field is missing.");
            ComputeInterfaceWeight(*WeightPairsA[j][k], *AvP_PairsA[j][k], *denomA, psi_PairsA[j][k].get());
        }
    }

    // -------------------------------------------------
    // Cathode particle-particle weights
    // -------------------------------------------------

    for (int j = 0;
         j < static_cast<int>(psC.size());
         ++j)
    {
        for (int k = j + 1;
             k < static_cast<int>(psC.size());
             ++k)
        {
            MFEM_VERIFY(WeightPairsC[j][k], "Cathode particle-pair weight field is missing.");
            ComputeInterfaceWeight(*WeightPairsC[j][k], *AvP_PairsC[j][k], *denomC, psi_PairsC[j][k].get());
        }
    }
}


void Domain_Parameters::CalculateTotals(const mfem::ParGridFunction &grid_function, const mfem::Vector &element_volumes, double &local_total, double &global_total)
{
    MFEM_VERIFY(element_volumes.Size() == pmesh->GetNE(), "Element-volume vector size does not match the local element count.");

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

        for (int j = 0;
             j < nodal_values.Size();
             ++j)
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
    MFEM_VERIFY(psi, "Half-cell total solid phase field is not initialized.");
    MFEM_VERIFY(ps.size() == particle_labels.size(), "Half-cell particle-field count does not match particle labels.");

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
    MFEM_VERIFY(psiA, "Full-cell anode phase field is not initialized.");
    MFEM_VERIFY(psiC, "Full-cell cathode phase field is not initialized.");
    MFEM_VERIFY(psA.size() == cfg.anode_materials.size(), "Anode material count does not match anode particle count.");
    MFEM_VERIFY(psC.size() == cfg.cathode_materials.size(), "Cathode material count does not match cathode particle count.");

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

    for (std::size_t k = 0;
         k < psC.size();
         ++k)
    {
        CalculateTotalPhaseField(*psC[k], tPsC[k], gtPsC[k]);
        CalculateTargetCurrent(tPsC[k], gTrgPsC[k], cfg.cathode_materials[k]);
        gTrgIC += gTrgPsC[k];
    }

    // gTrgIA = gTrgIC;
    // gTrgI = std::min(std::abs(gTrgIA), std::abs(gTrgIC));
    gTrgI = gTrgIC; // use cathode target current
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

