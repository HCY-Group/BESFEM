#include "../include/Initialize_Geometry.hpp"
#include "../include/Domain_Parameters.hpp"
#include "../include/MaterialProperties.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
void ClampPhase(mfem::ParGridFunction& phase)
{
    for (int i = 0; i < phase.Size(); ++i)
    {
        phase(i) = std::clamp(phase(i), 1.0e-6, 1.0);
    }
}

void PrintParticleTotals(const char* prefix, const std::vector<double>& totals,
    const std::vector<double>& currents)
{
    for (std::size_t k = 0; k < totals.size(); ++k)
    {
        std::cout << prefix << k << " phase total: " << totals[k]
            << ", target current: " << currents[k] << '\n';
    }
}

}

Domain_Parameters::Domain_Parameters(Initialize_Geometry &geo, const SimulationConfig &cfg)
    : geometry(geo), cfg(cfg),
    pmesh(geo.parallelMesh.get()), fespace(geo.parfespace),
    particle_labels(geo.particle_labels), anode_particle_labels(geo.anode_particle_labels), cathode_particle_labels(geo.cathode_particle_labels)
{}

Domain_Parameters::~Domain_Parameters() = default;

void Domain_Parameters::SetupDomainParameters()
{

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

void Domain_Parameters::InitializeGridFunctions()
{

    if (!fespace) {
        throw std::runtime_error("Finite element space is not initialized.");
    }

    for (auto* field : {&psi, &pse, &AvP, &AvB, &AvE, &denom})
    {
        *field = std::make_unique<mfem::ParGridFunction>(fespace.get());
    }

    if (cfg.mode == sim::CellMode::HALF)
    {
        InitializeHalfCellGridFunctions();
    }
    else
    {
        InitializeFullCellGridFunctions();
    }
}

Domain_Parameters::ParticleGroups Domain_Parameters::GetParticleGroups(sim::Electrode electrode)
{
    if (cfg.mode == sim::CellMode::HALF)
    {
        return {ps, AvPs, AvEs, WeightEs, AvP_Pairs, psi_Pairs, WeightPairs, tPs, gtPs, gTrgPs};
    }
    if (electrode == sim::Electrode::ANODE)
    {
        return {psA, AvPsA, AvEsA, WeightEsA, AvP_PairsA, psi_PairsA, WeightPairsA, tPsA, gtPsA, gTrgPsA};
    }
    MFEM_VERIFY(electrode == sim::Electrode::CATHODE, "Select one electrode's particle groups.");
    return {psC, AvPsC, AvEsC, WeightEsC, AvP_PairsC, psi_PairsC, WeightPairsC, tPsC, gtPsC, gTrgPsC};
}

void Domain_Parameters::AllocateParticleGroups(ParticleGroups groups, std::size_t count)
{
    for (auto* fields : {&groups.phase, &groups.gradient, &groups.electrolyte_interface, &groups.electrolyte_weight})
    {
        fields->clear();
        fields->resize(count);
        for (auto& field : *fields)
        {
            field = std::make_unique<mfem::ParGridFunction>(fespace.get());
        }
    }
    for (auto* pairs : {&groups.pair_interface, &groups.pair_phase, &groups.pair_weight})
    {
        pairs->clear();
        pairs->resize(count);
        for (std::size_t j = 0; j < count; ++j)
        {
            (*pairs)[j].resize(count);
            // All consumers address unordered pairs through j < k.
            for (std::size_t k = j + 1; k < count; ++k)
            {
                (*pairs)[j][k] = std::make_unique<mfem::ParGridFunction>(fespace.get());
            }
        }
    }
    groups.local_total.assign(count, 0.0);
    groups.global_total.assign(count, 0.0);
    groups.target_current.assign(count, 0.0);
}

void Domain_Parameters::CopyParticleMasks(ParticleGroups groups, const FieldList& masks,
    mfem::ParGridFunction& total)
{
    MFEM_VERIFY(masks.size() == groups.phase.size(), "One geometry mask is required per particle group.");
    total = 0.0;
    for (std::size_t k = 0; k < masks.size(); ++k)
    {
        *groups.phase[k] = *masks[k];
        total += *groups.phase[k];
        ClampPhase(*groups.phase[k]);
    }
}

void Domain_Parameters::BuildParticleInterfaces(ParticleGroups groups, mfem::ParGridFunction& denominator)
{
    const std::size_t count = groups.phase.size();
    for (std::size_t k = 0; k < count; ++k)
    {
        ComputeGradientMagnitude(*groups.phase[k], *groups.gradient[k]);
        BuildInterface(*groups.electrolyte_interface[k], *AvE, *groups.gradient[k]);
    }
    denominator = 0.0;
    for (std::size_t j = 0; j < count; ++j)
    {
        for (std::size_t k = j + 1; k < count; ++k)
        {
            BuildInterface(*groups.pair_interface[j][k], *groups.gradient[j], *groups.gradient[k]);
            BuildPairPhaseMask(*groups.pair_phase[j][k], *groups.phase[j], *groups.phase[k]);
            denominator += *groups.pair_interface[j][k];
        }
    }
    for (const auto& interface : groups.electrolyte_interface)
    {
        denominator += *interface;
    }
    for (std::size_t k = 0; k < count; ++k)
    {
        ComputeInterfaceWeight(*groups.electrolyte_weight[k], *groups.electrolyte_interface[k], denominator);
    }
    for (std::size_t j = 0; j < count; ++j)
    {
        for (std::size_t k = j + 1; k < count; ++k)
        {
            ComputeInterfaceWeight(*groups.pair_weight[j][k], *groups.pair_interface[j][k],
                denominator, groups.pair_phase[j][k].get());
        }
    }
}

double Domain_Parameters::CalculateParticleTotals(ParticleGroups groups,
    const std::vector<sim::MaterialType>& materials)
{
    MFEM_VERIFY(materials.size() == groups.phase.size(), "One material is required per particle group.");
    double target = 0.0;
    for (std::size_t k = 0; k < groups.phase.size(); ++k)
    {
        CalculateTotalPhaseField(*groups.phase[k], groups.local_total[k], groups.global_total[k]);
        CalculateTargetCurrent(groups.local_total[k], groups.target_current[k], materials[k]);
        target += groups.target_current[k];
    }
    return target;
}

void Domain_Parameters::InitializeHalfCellGridFunctions()
{
    AllocateParticleGroups(GetParticleGroups(cfg.half_electrode), particle_labels.size());
}

void Domain_Parameters::InitializeFullCellGridFunctions()
{
    for (auto* field : {&psiA, &psiC, &AvPA, &AvPC, &denomA, &denomC})
    {
        *field = std::make_unique<mfem::ParGridFunction>(fespace.get());
    }
    AllocateParticleGroups(GetParticleGroups(sim::Electrode::ANODE), anode_particle_labels.size());
    AllocateParticleGroups(GetParticleGroups(sim::Electrode::CATHODE), cathode_particle_labels.size());
}

void Domain_Parameters::InterpolateDomainParameters()
{

    if (cfg.mode == sim::CellMode::HALF)
    {
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
    *pse = *geometry.MaskFilterPse;
    CopyParticleMasks(GetParticleGroups(cfg.half_electrode), geometry.MaskFilters, *psi);
    ClampPhase(*psi);
    ClampPhase(*pse);
}

void Domain_Parameters::InterpolateFullCellMasks()
{
    *pse = *geometry.MaskFilterPse;
    CopyParticleMasks(GetParticleGroups(sim::Electrode::ANODE), geometry.MaskFiltersAnode, *psiA);
    CopyParticleMasks(GetParticleGroups(sim::Electrode::CATHODE), geometry.MaskFiltersCathode, *psiC);

    // Sum the raw electrode masks before clamping either electrode total.
    *psi = *psiA;
    *psi += *psiC;
    for (auto* phase : {psiA.get(), psiC.get(), psi.get(), pse.get()})
    {
        ClampPhase(*phase);
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

void Domain_Parameters::BuildInterface(mfem::ParGridFunction& out,
    const mfem::ParGridFunction& gradient_a, const mfem::ParGridFunction& gradient_b)
{
    out = gradient_a;
    out *= gradient_b;
    for (int i = 0; i < out.Size(); ++i)
    {
        out(i) = std::sqrt(out(i));
    }
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
    BuildParticleInterfaces(GetParticleGroups(cfg.half_electrode), *denom);

    for (std::size_t j = 0; j < ps.size(); ++j)
    {
        for (std::size_t k = j + 1; k < ps.size(); ++k)
        {
            const std::string filename = "AvP_Pair_" + std::to_string(j) + "_" + std::to_string(k);
            AvP_Pairs[j][k]->SaveAsOne(filename.c_str());
        }
    }
}

void Domain_Parameters::BuildFullCellInterfaces()
{
    ComputeGradientMagnitude(*psi, *AvP);
    ComputeGradientMagnitude(*psiA, *AvPA);
    ComputeGradientMagnitude(*psiC, *AvPC);
    ComputeGradientMagnitude(*pse, *AvE);
    BuildParticleInterfaces(GetParticleGroups(sim::Electrode::ANODE), *denomA);
    BuildParticleInterfaces(GetParticleGroups(sim::Electrode::CATHODE), *denomC);
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
    CalculateTotals(grid_function, EVol, local_total, global_total);
}

void Domain_Parameters::CalculatePhasePotentialsAndTargetCurrent()
{
    const int local_element_count = pmesh->GetNE();
    EVol.SetSize(local_element_count);

    for (int ei = 0; ei < local_element_count; ++ei)
    {
        EVol(ei) = pmesh->GetElementVolume(ei);
    }

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
    const auto& materials = cfg.half_electrode == sim::Electrode::CATHODE
        ? cfg.cathode_materials : cfg.anode_materials;
    CalculateTotalPhaseField(*psi, tPsi, gtPsi);
    gTrgI = CalculateParticleTotals(GetParticleGroups(cfg.half_electrode), materials);
}

void Domain_Parameters::CalculateFullCellPhasePotentialsAndTargetCurrent()
{
    CalculateTotalPhaseField(*psiA, tPsiA, gtPsiA);
    CalculateTotalPhaseField(*psiC, tPsiC, gtPsiC);
    CalculateTotalPhaseField(*psi, tPsi, gtPsi);
    gTrgIA = CalculateParticleTotals(GetParticleGroups(sim::Electrode::ANODE), cfg.anode_materials);
    gTrgIC = CalculateParticleTotals(GetParticleGroups(sim::Electrode::CATHODE), cfg.cathode_materials);
    gTrgI = gTrgIC; // Preserve the cathode-based target for full cells.
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

        PrintParticleTotals("Particle ", gtPs, gTrgPs);
    }
    else
    {
        std::cout << "Total anode phase: " << gtPsiA << '\n' << "Total cathode phase: " << gtPsiC << '\n'
            << "Anode capacity-based current: " << gTrgIA << '\n' << "Cathode capacity-based current: " << gTrgIC << '\n'
            << "Selected full-cell target current: " << gTrgI << '\n';

        PrintParticleTotals("Anode particle ", gtPsA, gTrgPsA);

        PrintParticleTotals("Cathode particle ", gtPsC, gTrgPsC);
    }
}

