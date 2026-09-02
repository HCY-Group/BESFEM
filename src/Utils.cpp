#include "../include/Utils.hpp"
#include "../include/Constants.hpp"
#include "../include/SimulationConfig.hpp"
#include "../include/MaterialProperties.hpp"
#include "../include/SimulationState.hpp"

#include <numeric>

namespace
{

    template <typename ParticleState>
    void BuildCombinedParticleConcentration(const std::vector<ParticleState>& particles, const std::vector<std::unique_ptr<mfem::ParGridFunction>>& particle_ps, mfem::ParGridFunction& output)
    {
        const int np = static_cast<int>(particles.size());
        output = 0.0;
        mfem::ParGridFunction denominator(output.ParFESpace());
        denominator = 0.0;

        for (int k = 0; k < np; ++k)
        {
            for (int i = 0; i < output.Size(); ++i)
            {
                const double psi = (*particle_ps[k])(i);
                output(i) += psi * (*particles[k].Cn_gf)(i);
                denominator(i) += psi;
            }
        }

        constexpr double eps = 1.0e-30;

        for (int i = 0; i < output.Size(); ++i)
        {
            const double d = denominator(i);

            if (d > eps)
            {
                output(i) /= d;
            }
            else
            {
                output(i) = 0.0;
            }
        }
    }

    std::string MakeTimestepSuffix(int t)
    {
        std::ostringstream step;
        step << "_" << std::setw(7) << std::setfill('0') << t;
        return step.str();
    }
} 

Utils::Utils(Initialize_Geometry &geo, Domain_Parameters &para, const SimulationConfig &cfg)
    : geometry_(geo), domain_(para), cfg(cfg), pmesh_(geo.parallelMesh.get()), fes_(geo.parfespace),
      nE_(geo.nE), nC_(geo.nC), nV_(geo.nV), EVol_(para.EVol), EAvg_(geo.nE), VtxVal_(geo.nC), TmpF_(geo.parfespace.get())
{
}

void Utils::InitializeReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value)
{
    Rx2 = Rx1;
    Rx2 *= value;
}

void Utils::InitializeReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &Rx3, double value)
{
    Rx3 = Rx1;
    Rx3 += Rx2;
    Rx3 *= value;
}

void Utils::CalculateLithiation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, double gtps)
{
    TmpF_ = Cn;
    TmpF_ *= psx;

    double local_sum = 0.0;

    for (int ei = 0; ei < nE_; ei++)
    {
        TmpF_.GetNodalValues(ei, VtxVal_);
        double sum = std::accumulate(VtxVal_.begin(), VtxVal_.end(), 0.0);

        EAvg_(ei) = sum / nC_;
        local_sum += EAvg_(ei) * EVol_(ei);
    }

    MPI_Allreduce(&local_sum, &Xfr_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Xfr_ /= gtps;
}

void Utils::CalculateReactionInfx(mfem::ParGridFunction &Rx, double &infx_)
{
    double xCrnt = 0.0;

    mfem::Vector Rmin, Rmax;
    pmesh_->GetBoundingBox(Rmin, Rmax);

    double Lw = (Rmax(1) - Rmin(1));
    if (pmesh_->Dimension() == 3)
        Lw *= (Rmax(2) - Rmin(2));

    for (int ei = 0; ei < nE_; ei++)
    {
        Rx.GetNodalValues(ei, VtxVal_);
        double sum = 0;
        for (int i = 0; i < nC_; i++) sum += VtxVal_[i];

        EAvg_(ei) = sum / nC_;
        xCrnt += EAvg_(ei) * EVol_(ei);
    }

    MPI_Allreduce(&xCrnt, &geCrnt_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    infx_ = geCrnt_ / Lw;
}

void Utils::CalculateGlobalError(mfem::ParGridFunction &px0, mfem::ParGridFunction &potential, mfem::ParGridFunction &psx, double &globalerror, double gtPsx)
{
    TmpF_ = px0;
    TmpF_ -= potential;
    TmpF_ *= TmpF_;
    TmpF_ *= psx;

    double local_sum = TmpF_.Sum();
    MPI_Allreduce(&local_sum, &globalerror, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    globalerror /= gtPsx;
}

void Utils::ComputePairFlux(mfem::ParGridFunction &sum_part, mfem::ParGridFunction &weight, mfem::ParGridFunction &grad_psi, mfem::ParGridFunction &mu_1, mfem::ParGridFunction &mu_2, double rho)
{
    for (int vi = 0; vi < nV_; vi++){

        double grad_psi_val = grad_psi(vi);
        double weight_val = weight(vi);
        double mu1_val = mu_1(vi);
        double mu2_val = mu_2(vi);

        // const double rho = MaterialProperties::SiteDensity(cfg.cathode_materials[0]);

        sum_part(vi) = weight_val * grad_psi_val * rho * (1.0/Constants::RT) * Constants::Perm * (mu2_val - mu1_val);
    }

}

void Utils::SetInitialValue(mfem::ParGridFunction &Cn, double initial_value)
    {
        for (int i = 0; i < Cn.Size(); i++)
            Cn(i) = initial_value;
    }


void Utils::SaveHalfCellSnapshot(int t, const std::string& outdir, Initialize_Geometry& geometry,
    Domain_Parameters& domain_parameters, SimulationState& state, sim::Electrode electrode, int save_interval)
{

    if (t % save_interval != 0)
    {
        return;
    }

    const std::string suffix = MakeTimestepSuffix(t);

    if (t == 0)
    {
        geometry.parallelMesh->SaveAsOne((outdir + "/pmesh").c_str());
        domain_parameters.psi->SaveAsOne((outdir + "/psi").c_str());
        domain_parameters.pse->SaveAsOne((outdir + "/pse").c_str());

        // Save each individual particle/material phase field.
        for (int k = 0; k < static_cast<int>(domain_parameters.ps.size()); ++k)
        {
            std::ostringstream name;
            name << outdir << "/ps_" << k;
            domain_parameters.ps[k]->SaveAsOne(name.str().c_str());
        }
    }

    state.CnE_gf->SaveAsOne((outdir + "/CnE" + suffix).c_str());
    state.phE_gf->SaveAsOne((outdir + "/phE" + suffix).c_str());

    mfem::ParGridFunction CnP(geometry.parfespace.get());

    if (electrode == sim::Electrode::ANODE)
    {
        BuildCombinedParticleConcentration(state.anode_particles, domain_parameters.ps, CnP);
        CnP.SaveAsOne((outdir + "/CnP" + suffix).c_str());
        state.phA_gf->SaveAsOne((outdir + "/phP" + suffix).c_str());
    }
    else
    {
        BuildCombinedParticleConcentration(state.cathode_particles, domain_parameters.ps, CnP);
        CnP.SaveAsOne((outdir + "/CnP" + suffix).c_str());
        state.phC_gf->SaveAsOne((outdir + "/phP" + suffix).c_str());
    }

    state.Rxn_gf->SaveAsOne((outdir + "/RxnP" + suffix).c_str());
}

void Utils::SaveFullCellSnapshot(int t, const std::string& outdir, Initialize_Geometry& geometry,
    Domain_Parameters& domain_parameters, SimulationState& state, int save_interval)
{
    if (t % save_interval != 0)
    {
        return;
    }

    const std::string suffix = MakeTimestepSuffix(t);

    if (t == 0)
    {
        geometry.parallelMesh->SaveAsOne((outdir + "/pmesh").c_str());
        domain_parameters.psiA->SaveAsOne((outdir + "/psiA").c_str());
        domain_parameters.psiC->SaveAsOne((outdir + "/psiC").c_str());
        domain_parameters.pse->SaveAsOne((outdir + "/pse").c_str());

        for (int k = 0; k < static_cast<int>(domain_parameters.psA.size()); ++k)
        {
            std::ostringstream name;
            name << outdir << "/psA_" << k;
            domain_parameters.psA[k]->SaveAsOne(name.str().c_str());
        }

        for (int k = 0; k < static_cast<int>(domain_parameters.psC.size()); ++k)
        {

            std::ostringstream name;
            name << outdir << "/psC_" << k;
            domain_parameters.psC[k]->SaveAsOne(name.str().c_str());
        }
    }

    mfem::ParGridFunction CnA(geometry.parfespace.get());
    mfem::ParGridFunction CnC(geometry.parfespace.get());

    BuildCombinedParticleConcentration(state.anode_particles, domain_parameters.psA, CnA);
    BuildCombinedParticleConcentration(state.cathode_particles, domain_parameters.psC, CnC);

    CnA.SaveAsOne((outdir + "/CnA" + suffix).c_str());
    CnC.SaveAsOne((outdir + "/CnC" + suffix).c_str());
    state.CnE_gf->SaveAsOne((outdir + "/CnE" + suffix).c_str());

    state.phA_gf->SaveAsOne((outdir + "/phA" + suffix).c_str());
    state.phC_gf->SaveAsOne((outdir + "/phC" + suffix).c_str());
    state.phE_gf->SaveAsOne((outdir + "/phE" + suffix).c_str());

    state.RxnA_gf->SaveAsOne((outdir + "/RxnA" + suffix).c_str());
    state.RxnC_gf->SaveAsOne((outdir + "/RxnC" + suffix).c_str());
    state.RxnE_gf->SaveAsOne((outdir + "/RxnE" + suffix).c_str());
}

void Utils::PrintSimulationParameters(const SimulationConfig &cfg, const std::string &outdir)
{
    if (mfem::Mpi::WorldRank() != 0)
    {
        return;
    }

    std::cout << "\n===== Simulation Parameters =====\n"
              << "output_dir = " << outdir << "\n"
              << "dt   = " << cfg.dt << "\n"
              << "dh   = " << cfg.dh << "\n"
              << "gc   = " << cfg.gc << "\n"
              << "Cr   = " << cfg.Cr << "\n"
              << "Vsr0 = " << cfg.Vsr0 << "\n"
              << "=================================\n"
              << std::endl;
}

void Utils::PrintHalfCellStatus(int t, double VCell,double total_current, double total_target,
    const std::vector<double> &particle_currents, const SimulationState &state, const Domain_Parameters &para, sim::Electrode electrode)
{
    if (mfem::Mpi::WorldRank() != 0)
    {
        return;
    }

    // const auto &particles = (electrode == sim::Electrode::ANODE) ? state.anode_particles : state.cathode_particles;
    const bool is_anode = (electrode == sim::Electrode::ANODE);
    const int np = is_anode ? state.anode_particles.size() : state.cathode_particles.size();

    std::cout << "timestep: " << t << ", VCell = " << VCell << ", TotalCurrent = " << total_current << ", TotalTarget = " << total_target;

    for (int j = 0; j < np; ++j)
    {
        std::cout << ", Current_" << j << " = " << particle_currents[j] << ", Target_" << j << " = " << para.gTrgPs[j];
    }

    std::cout << std::endl;

    const char *electrode_name = (electrode == sim::Electrode::ANODE) ? "ANODE" : "CATHODE";

    std::cout << "timestep: " << t << " [" << electrode_name << " HALF-CELL]" << ", VCell = " << VCell << ", BvE = " << state.electrolyte_potential->GetBoundaryVoltage();

    if (np > 0)
    {
        if (is_anode)
        {
            std::cout << ", Cp_min = " << state.anode_particles[0].Cn_gf->Min()
                    << ", Cp_max = " << state.anode_particles[0].Cn_gf->Max();
        }
        else
        {
            std::cout << ", Cp_min = " << state.cathode_particles[0].Cn_gf->Min()
                    << ", Cp_max = " << state.cathode_particles[0].Cn_gf->Max();
        }
    }

    std::cout << ", Ce_min = " << state.CnE_gf->Min() << ", Ce_max = " << state.CnE_gf->Max() << std::endl;

    double Xfr_avg = 0.0;
    double total_weight = 0.0;

    for (int j = 0; j < np; ++j)
    {

        double Xfr_j;
        if (is_anode)
        {
            Xfr_j = state.anode_particles[j].concentration->GetLithiation();
        }
        else
        {
            Xfr_j = state.cathode_particles[j].concentration->GetLithiation();
        }
        const double weight_j = para.gtPs[j];

        Xfr_avg += weight_j * Xfr_j;
        total_weight += weight_j;

        std::cout << ", Xfr_" << j << " = " << Xfr_j;
    }

    if (total_weight > 0.0)
    {
        Xfr_avg /= total_weight;
    }

    std::cout << ", Xfr"  << (electrode == sim::Electrode::ANODE ? "A" : "C") << "_avg = " << Xfr_avg << std::endl;

    const double XfrE = state.electrolyte_concentration->GetLithiation();
    const double weight_E = para.gtPse;
    const double XfrE_avg = (weight_E > 0.0) ? XfrE : 0.0;

    std::cout << "XfrE = " << XfrE << ", weight_je = " << weight_E << ", XfrE_avg = " << XfrE_avg << std::endl;
}

void Utils::PrintFullCellStatus(int t, double VCell, double anode_current, double cathode_current, const SimulationState &state, const Domain_Parameters &para)
{
    if (mfem::Mpi::WorldRank() != 0)
    {
        return;
    }

    const int npA = static_cast<int>(state.anode_particles.size());
    const int npC = static_cast<int>(state.cathode_particles.size());

    double XfrA_avg = 0.0;
    double total_anode_weight = 0.0;

    for (int j = 0; j < npA; ++j)
    {
        const double Xfr_j = state.anode_particles[j].concentration->GetLithiation();
        const double weight_j = para.gtPsA[j];

        XfrA_avg += weight_j * Xfr_j;
        total_anode_weight += weight_j;
    }

    if (total_anode_weight > 0.0)
    {
        XfrA_avg /= total_anode_weight;
    }

    // ------------------------------------------------------------
    // Calculate average cathode lithiation
    // ------------------------------------------------------------

    double XfrC_avg = 0.0;
    double total_cathode_weight = 0.0;

    for (int j = 0; j < npC; ++j)
    {
        const double Xfr_j = state.cathode_particles[j].concentration->GetLithiation();
        const double weight_j = para.gtPsC[j];

        XfrC_avg += weight_j * Xfr_j;
        total_cathode_weight += weight_j;
    }

    if (total_cathode_weight > 0.0)
    {
        XfrC_avg /= total_cathode_weight;
    }

    std::cout << "timestep: " << t << " [FULL-CELL]"
              << ", XfrA_avg = " << XfrA_avg << ", XfrC_avg = " << XfrC_avg
              << ", Anode current = " << anode_current << ", Cathode current = " << cathode_current
              << ", Anode target = " << para.gTrgI  << ", Cathode target = " << para.gTrgI
              << ", VCell = " << VCell << ", BvA = "  << state.anode_potential->GetBoundaryVoltage()
              << ", BvC = " << state.cathode_potential->GetBoundaryVoltage()
              << ", BvE = " << state.electrolyte_potential->GetBoundaryVoltage() << std::endl;


    for (int j = 0; j < npA; ++j)
    {
        std::cout << "    Anode particle " << j << ", Xfr = " << state.anode_particles[j].concentration->GetLithiation() << std::endl;
    }

    for (int j = 0; j < npC; ++j)
    {
        std::cout << "    Cathode particle " << j << ", Xfr = " << state.cathode_particles[j].concentration->GetLithiation() << std::endl;
    }
}

void Utils::PrintProgramTime(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
    if (mfem::Mpi::WorldRank() != 0)
    {
        return;
    }

    const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "Total Program Time: " << elapsed.count() << " seconds" << std::endl;
}