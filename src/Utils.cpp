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

        MFEM_VERIFY(static_cast<int>(particle_ps.size()) == np, "BuildCombinedParticleConcentration: particle and phase-field counts differ.");

        output = 0.0;

        mfem::ParGridFunction denominator(output.ParFESpace());
        denominator = 0.0;

        for (int k = 0; k < np; ++k)
        {
            MFEM_VERIFY(particles[k].Cn_gf, "BuildCombinedParticleConcentration: concentration field is null.");
            MFEM_VERIFY(particle_ps[k], "BuildCombinedParticleConcentration: phase field is null.");
            MFEM_VERIFY(particles[k].Cn_gf->Size() == output.Size(), "BuildCombinedParticleConcentration: concentration field size mismatch.");
            MFEM_VERIFY(particle_ps[k]->Size() == output.Size(), "BuildCombinedParticleConcentration: phase field size mismatch.");

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

// Full Cell
// void Utils::SaveSimulationSnapshot(int t, const std::string &outdir, Initialize_Geometry &geometry, Domain_Parameters &domain_parameters, mfem::ParGridFunction &phA,
//     mfem::ParGridFunction &phC, mfem::ParGridFunction &phE, mfem::ParGridFunction &CnA, mfem::ParGridFunction &CnC, mfem::ParGridFunction &CnE, mfem::ParGridFunction &CnApsi,
//     mfem::ParGridFunction &CnCpsi, mfem::ParGridFunction &CnEpsi, mfem::ParGridFunction &CnP, int save_interval)
// {
//     if (t % save_interval != 0) return;

//     std::ostringstream step;
//     step << "_" << std::setw(5) << std::setfill('0') << t;
//     std::string suff = step.str();

//     if (t == 0)
//     {
//         geometry.parallelMesh->SaveAsOne((outdir + "/pmesh").c_str());
//         domain_parameters.psiA->SaveAsOne((outdir + "/psiA").c_str());
//         domain_parameters.psiC->SaveAsOne((outdir + "/psiC").c_str());
//         domain_parameters.pse->SaveAsOne((outdir + "/pse").c_str());
//     }

//     phA.SaveAsOne((outdir + "/phA" + suff).c_str());
//     phC.SaveAsOne((outdir + "/phC" + suff).c_str());
//     phE.SaveAsOne((outdir + "/phE" + suff).c_str());

//     CnA.SaveAsOne((outdir + "/CnA_raw" + suff).c_str());
//     CnC.SaveAsOne((outdir + "/CnC_raw" + suff).c_str());
//     CnE.SaveAsOne((outdir + "/CnE_raw" + suff).c_str());

//     CnApsi = CnA; CnApsi *= *domain_parameters.psiA;
//     CnApsi.SaveAsOne((outdir + "/CnA" + suff).c_str());

//     CnCpsi = CnC; CnCpsi *= *domain_parameters.psiC;
//     CnCpsi.SaveAsOne((outdir + "/CnC" + suff).c_str());

//     CnEpsi = CnE; CnEpsi *= *domain_parameters.pse;
//     CnEpsi.SaveAsOne((outdir + "/CnE" + suff).c_str());

//     CnP = CnApsi;
//     CnP += CnCpsi;
//     CnP.SaveAsOne((outdir + "/CnP" + suff).c_str());
// }

// // Half Cell
// void Utils::SaveSimulationSnapshot(int t, const std::string &outdir,
//     Initialize_Geometry &geometry, Domain_Parameters &domain_parameters, mfem::ParGridFunction &phC, mfem::ParGridFunction &phE, mfem::ParGridFunction &CnC, mfem::ParGridFunction &CnE,
//     mfem::ParGridFunction &CnCpsi, mfem::ParGridFunction &CnEpsi, int save_interval)
// {
//     if (t % save_interval != 0) return;

//     std::ostringstream step;
//     step << "_" << std::setw(5) << std::setfill('0') << t;
//     std::string suff = step.str();

//     if (t == 0)
//     {
//         geometry.parallelMesh->SaveAsOne((outdir + "/pmesh").c_str());
//         domain_parameters.psi->SaveAsOne((outdir + "/psi").c_str());
//         domain_parameters.pse->SaveAsOne((outdir + "/pse").c_str());
//     }

//     phC.SaveAsOne((outdir + "/phC" + suff).c_str());
//     phE.SaveAsOne((outdir + "/phE" + suff).c_str());

//     CnC.SaveAsOne((outdir + "/CnP_raw" + suff).c_str());
//     CnE.SaveAsOne((outdir + "/CnE_raw" + suff).c_str());

//     CnCpsi = CnC; CnCpsi *= *domain_parameters.psi;
//     CnCpsi.SaveAsOne((outdir + "/CnP" + suff).c_str());

//     CnEpsi = CnE; CnEpsi *= *domain_parameters.pse;
//     CnEpsi.SaveAsOne((outdir + "/CnE" + suff).c_str());
// }

void Utils::SetInitialValue(mfem::ParGridFunction &Cn, double initial_value)
    {
        for (int i = 0; i < Cn.Size(); i++)
            Cn(i) = initial_value;
    }


// void Utils::SaveSimulationSnapshotMulti(int t, const std::string& outdir, Initialize_Geometry& geometry,
//     Domain_Parameters& domain_parameters, const std::vector<mfem::ParGridFunction*>& particle_cn,
//     const std::vector<std::unique_ptr<mfem::ParGridFunction>>& particle_ps, mfem::ParGridFunction& electrode_psi,
//     std::vector<std::unique_ptr<mfem::ParGridFunction>>& particle_out, const std::string& electrode_name, int save_interval)
//     {
//         if (t % save_interval != 0)
//         {
//             return;
//         }

//         const int np = static_cast<int>(particle_cn.size());
//         MFEM_VERIFY(static_cast<int>(particle_ps.size()) == np, "SaveSimulationSnapshotMulti: particle concentration and phase-field counts differ.");
//         MFEM_VERIFY(static_cast<int>(particle_out.size()) == np, "SaveSimulationSnapshotMulti: particle concentration and output counts differ.");

//         std::ostringstream step;
//         step << "_" << std::setw(5) << std::setfill('0') << t;
//         const std::string suff = step.str();

//         // =====================================================================
//         // Save mesh and phase fields at timestep zero
//         // =====================================================================

//         if (t == 0)
//         {
//             MFEM_VERIFY(geometry.parallelMesh, "SaveSimulationSnapshotMulti: parallel mesh is null.");
//             MFEM_VERIFY(domain_parameters.pse, "SaveSimulationSnapshotMulti: electrolyte phase field is null.");

//             geometry.parallelMesh->SaveAsOne( (outdir + "/pmesh").c_str());
//             electrode_psi.SaveAsOne((outdir + "/psi" + electrode_name).c_str());
//             domain_parameters.pse->SaveAsOne((outdir + "/pse").c_str());

//             for (int k = 0; k < np; ++k)
//             {
//                 MFEM_VERIFY(particle_ps[k], "SaveSimulationSnapshotMulti: null particle phase field.");
//                 std::ostringstream psi_name;
//                 psi_name << outdir << "/ps" << electrode_name << "_" << (k + 1);
//                 particle_ps[k]->SaveAsOne(psi_name.str().c_str());
//             }
//         }

//         // =====================================================================
//         // Save each particle concentration and masked concentration
//         // =====================================================================

//         for (int k = 0; k < np; ++k)
//         {
//             MFEM_VERIFY(particle_cn[k] != nullptr, "SaveSimulationSnapshotMulti: null particle concentration field.");
//             MFEM_VERIFY(particle_ps[k], "SaveSimulationSnapshotMulti: null particle phase field.");
//             MFEM_VERIFY(particle_out[k], "SaveSimulationSnapshotMulti: null particle output field.");
//             MFEM_VERIFY(particle_cn[k]->Size() == particle_ps[k]->Size(), "SaveSimulationSnapshotMulti: concentration and phase-field sizes differ.");
//             MFEM_VERIFY(particle_out[k]->Size() == particle_cn[k]->Size(), "SaveSimulationSnapshotMulti: output and concentration sizes differ.");

//             std::ostringstream raw_name;
//             std::ostringstream masked_name;

//             // raw_name << outdir << "/Cn" << electrode_name << "_" << (k + 1) << suff;
//             // masked_name << outdir << "/C" << electrode_name << "_" << (k + 1) << "_out" << suff;

//             // particle_cn[k]->SaveAsOne(raw_name.str().c_str());
//             *particle_out[k] = *particle_cn[k];
//             *particle_out[k] *= *particle_ps[k];
//             // particle_out[k]->SaveAsOne(masked_name.str().c_str());
//         }

//         // =====================================================================
//         // Build union mask and total particle concentration
//         // =====================================================================

//         mfem::ParGridFunction psi_union(geometry.parfespace.get());
//         mfem::ParGridFunction denom(geometry.parfespace.get());
//         mfem::ParGridFunction CnP_total(geometry.parfespace.get());

//         psi_union = 0.0;
//         denom = 0.0;
//         CnP_total = 0.0;

//         for (int k = 0; k < np; ++k)
//         {
//             psi_union += *particle_ps[k];
//             denom += *particle_ps[k];
//         }

//         for (int i = 0; i < psi_union.Size(); ++i)
//         {
//             psi_union(i) = std::min(1.0, psi_union(i));
//         }

//         const double eps = 1.0e-30;

//         for (int i = 0; i < denom.Size(); ++i)
//         {
//             const double d = denom(i);

//             if (d > eps)
//             {
//                 double numerator = 0.0;

//                 for (int k = 0; k < np; ++k)
//                 {
//                     numerator += (*particle_ps[k])(i) * (*particle_cn[k])(i);
//                 }

//                 CnP_total(i) = numerator / (d + eps);
//             }
//             else
//             {
//                 CnP_total(i) = 0.0;
//             }
//         }

//         CnP_total *= psi_union;
//         CnP_total.SaveAsOne((outdir + "/Cn" + electrode_name + "_total" + suff).c_str());
//     }


// void Utils::SaveCombinedSnapshot(
//     int t,
//     const std::string& outdir,
//     Initialize_Geometry& geometry,
//     Domain_Parameters& domain_parameters,
//     const std::vector<mfem::ParGridFunction*>& anode_fields,
//     const std::vector<mfem::ParGridFunction*>& anode_ps,
//     const std::vector<mfem::ParGridFunction*>& cathode_fields,
//     const std::vector<mfem::ParGridFunction*>& cathode_ps,
//     const std::vector<mfem::ParGridFunction*>& electrolyte_fields,
//     const std::vector<mfem::ParGridFunction*>& electrolyte_ps,
//     const std::string& filename,
//     int save_interval)
// {
//     if (t % save_interval != 0)
//     {
//         return;
//     }

//     if (t == 0)
//     {
//         MFEM_VERIFY(geometry.parallelMesh, "SaveSimulationSnapshotMulti: parallel mesh is null.");
//         MFEM_VERIFY(domain_parameters.pse, "SaveSimulationSnapshotMulti: electrolyte phase field is null.");

//         geometry.parallelMesh->SaveAsOne( (outdir + "/pmesh").c_str());
//         domain_parameters.psi->SaveAsOne((outdir + "/psi").c_str());
//         domain_parameters.pse->SaveAsOne((outdir + "/pse").c_str());
//     }

//     MFEM_VERIFY(
//         anode_fields.size() == anode_ps.size(),
//         "SaveCombinedSnapshot: anode field and phase-field counts differ.");

//     MFEM_VERIFY(
//         cathode_fields.size() == cathode_ps.size(),
//         "SaveCombinedSnapshot: cathode field and phase-field counts differ.");

//     MFEM_VERIFY(
//         electrolyte_fields.size() == electrolyte_ps.size(),
//         "SaveCombinedSnapshot: electrolyte field and phase-field counts differ.");

//     MFEM_VERIFY(
//         geometry.parfespace,
//         "SaveCombinedSnapshot: finite-element space is null.");

//     // ---------------------------------------------------------------------
//     // Time suffix
//     // ---------------------------------------------------------------------

//     std::ostringstream step;
//     step << "_" << std::setw(5) << std::setfill('0') << t;
//     const std::string suffix = step.str();

//     // ---------------------------------------------------------------------
//     // Combined quantities
//     // ---------------------------------------------------------------------

//     mfem::ParGridFunction numerator(geometry.parfespace.get());
//     mfem::ParGridFunction denominator(geometry.parfespace.get());
//     mfem::ParGridFunction combined_field(geometry.parfespace.get());

//     numerator = 0.0;
//     denominator = 0.0;
//     combined_field = 0.0;

//     // =====================================================================
//     // ANODE
//     // =====================================================================

//     for (std::size_t k = 0; k < anode_fields.size(); ++k)
//     {
//         MFEM_VERIFY(
//             anode_fields[k] != nullptr,
//             "SaveCombinedSnapshot: null anode field.");

//         MFEM_VERIFY(
//             anode_ps[k] != nullptr,
//             "SaveCombinedSnapshot: null anode phase field.");

//         MFEM_VERIFY(
//             anode_fields[k]->Size() == anode_ps[k]->Size(),
//             "SaveCombinedSnapshot: anode field and phase-field sizes differ.");

//         for (int i = 0; i < numerator.Size(); ++i)
//         {
//             const double psi = (*anode_ps[k])(i);

//             numerator(i) += psi * (*anode_fields[k])(i);
//             denominator(i) += psi;
//         }
//     }

//     // =====================================================================
//     // CATHODE
//     // =====================================================================

//     for (std::size_t k = 0; k < cathode_fields.size(); ++k)
//     {
//         MFEM_VERIFY(
//             cathode_fields[k] != nullptr,
//             "SaveCombinedSnapshot: null cathode field.");

//         MFEM_VERIFY(
//             cathode_ps[k] != nullptr,
//             "SaveCombinedSnapshot: null cathode phase field.");

//         MFEM_VERIFY(
//             cathode_fields[k]->Size() == cathode_ps[k]->Size(),
//             "SaveCombinedSnapshot: cathode field and phase-field sizes differ.");

//         for (int i = 0; i < numerator.Size(); ++i)
//         {
//             const double psi = (*cathode_ps[k])(i);

//             numerator(i) += psi * (*cathode_fields[k])(i);
//             denominator(i) += psi;
//         }
//     }

//     // =====================================================================
//     // ELECTROLYTE
//     // =====================================================================

//     for (std::size_t k = 0; k < electrolyte_fields.size(); ++k)
//     {
//         MFEM_VERIFY(
//             electrolyte_fields[k] != nullptr,
//             "SaveCombinedSnapshot: null electrolyte field.");

//         MFEM_VERIFY(
//             electrolyte_ps[k] != nullptr,
//             "SaveCombinedSnapshot: null electrolyte phase field.");

//         MFEM_VERIFY(
//             electrolyte_fields[k]->Size() == electrolyte_ps[k]->Size(),
//             "SaveCombinedSnapshot: electrolyte field and phase-field sizes differ.");

//         for (int i = 0; i < numerator.Size(); ++i)
//         {
//             const double psi = (*electrolyte_ps[k])(i);

//             numerator(i) += psi * (*electrolyte_fields[k])(i);
//             denominator(i) += psi;
//         }
//     }

//     // =====================================================================
//     // Construct final combined field
//     // =====================================================================

//     const double eps = 1.0e-30;

//     for (int i = 0; i < combined_field.Size(); ++i)
//     {
//         if (denominator(i) > eps)
//         {
//             combined_field(i) = numerator(i) / denominator(i);
//         }
//         else
//         {
//             combined_field(i) = 0.0;
//         }
//     }

//     // =====================================================================
//     // Save
//     // =====================================================================

//     combined_field.SaveAsOne(
//         (outdir + "/" + filename + suffix).c_str());
// }

void Utils::SaveHalfCellSnapshot(
    int t,
    const std::string& outdir,
    Initialize_Geometry& geometry,
    Domain_Parameters& domain_parameters,
    SimulationState& state,
    sim::Electrode electrode,
    int save_interval)
{
    // ========================================================================
    // SAVE INTERVAL
    // ========================================================================

    if (t % save_interval != 0)
    {
        return;
    }

    MFEM_VERIFY(
        geometry.parfespace,
        "SaveHalfCellSnapshot: finite element space is null.");

    MFEM_VERIFY(
        geometry.parallelMesh,
        "SaveHalfCellSnapshot: parallel mesh is null.");

    MFEM_VERIFY(
        domain_parameters.psi,
        "SaveHalfCellSnapshot: electrode phase field is null.");

    MFEM_VERIFY(
        domain_parameters.pse,
        "SaveHalfCellSnapshot: electrolyte phase field is null.");

    MFEM_VERIFY(
        state.CnE_gf,
        "SaveHalfCellSnapshot: electrolyte concentration is null.");

    MFEM_VERIFY(
        state.phE_gf,
        "SaveHalfCellSnapshot: electrolyte potential is null.");

    MFEM_VERIFY(
        state.Rxn_gf,
        "SaveHalfCellSnapshot: reaction field is null.");

    const std::string suffix = MakeTimestepSuffix(t);


    // ========================================================================
    // SAVE STATIC GEOMETRY AT t = 0
    // ========================================================================

    if (t == 0)
    {
        geometry.parallelMesh->SaveAsOne(
            (outdir + "/pmesh").c_str());

        domain_parameters.psi->SaveAsOne(
            (outdir + "/psi").c_str());

        domain_parameters.pse->SaveAsOne(
            (outdir + "/pse").c_str());

        // Save each individual particle/material phase field.
        for (int k = 0;
             k < static_cast<int>(domain_parameters.ps.size());
             ++k)
        {
            MFEM_VERIFY(
                domain_parameters.ps[k],
                "SaveHalfCellSnapshot: particle phase field is null.");

            std::ostringstream name;

            name << outdir
                 << "/ps_"
                 << k;

            domain_parameters.ps[k]->SaveAsOne(
                name.str().c_str());
        }
    }


    // ========================================================================
    // ELECTROLYTE
    // ========================================================================

    state.CnE_gf->SaveAsOne(
        (outdir + "/CnE" + suffix).c_str());

    state.phE_gf->SaveAsOne(
        (outdir + "/phE" + suffix).c_str());


    // ========================================================================
    // ACTIVE ELECTRODE
    // ========================================================================

    mfem::ParGridFunction CnP(
        geometry.parfespace.get());

    if (electrode == sim::Electrode::ANODE)
    {
        MFEM_VERIFY(
            state.phA_gf,
            "SaveHalfCellSnapshot: anode potential is null.");

        MFEM_VERIFY(
            !state.anode_particles.empty(),
            "SaveHalfCellSnapshot: no anode particles exist.");

        BuildCombinedParticleConcentration(
            state.anode_particles,
            domain_parameters.ps,
            CnP);

        CnP.SaveAsOne(
            (outdir + "/CnP" + suffix).c_str());

        state.phA_gf->SaveAsOne(
            (outdir + "/phP" + suffix).c_str());
    }
    else
    {
        MFEM_VERIFY(
            state.phC_gf,
            "SaveHalfCellSnapshot: cathode potential is null.");

        MFEM_VERIFY(
            !state.cathode_particles.empty(),
            "SaveHalfCellSnapshot: no cathode particles exist.");

        BuildCombinedParticleConcentration(
            state.cathode_particles,
            domain_parameters.ps,
            CnP);

        CnP.SaveAsOne(
            (outdir + "/CnP" + suffix).c_str());

        state.phC_gf->SaveAsOne(
            (outdir + "/phP" + suffix).c_str());
    }


    // ========================================================================
    // REACTION
    // ========================================================================

    state.Rxn_gf->SaveAsOne(
        (outdir + "/RxnP" + suffix).c_str());
}

void Utils::SaveFullCellSnapshot(
    int t,
    const std::string& outdir,
    Initialize_Geometry& geometry,
    Domain_Parameters& domain_parameters,
    SimulationState& state,
    int save_interval)
{
    // ========================================================================
    // SAVE INTERVAL
    // ========================================================================

    if (t % save_interval != 0)
    {
        return;
    }

    MFEM_VERIFY(
        geometry.parfespace,
        "SaveFullCellSnapshot: finite element space is null.");

    MFEM_VERIFY(
        geometry.parallelMesh,
        "SaveFullCellSnapshot: parallel mesh is null.");

    MFEM_VERIFY(
        domain_parameters.psiA,
        "SaveFullCellSnapshot: anode phase field is null.");

    MFEM_VERIFY(
        domain_parameters.psiC,
        "SaveFullCellSnapshot: cathode phase field is null.");

    MFEM_VERIFY(
        domain_parameters.pse,
        "SaveFullCellSnapshot: electrolyte phase field is null.");

    MFEM_VERIFY(
        state.CnE_gf,
        "SaveFullCellSnapshot: electrolyte concentration is null.");

    MFEM_VERIFY(
        state.phA_gf,
        "SaveFullCellSnapshot: anode potential is null.");

    MFEM_VERIFY(
        state.phC_gf,
        "SaveFullCellSnapshot: cathode potential is null.");

    MFEM_VERIFY(
        state.phE_gf,
        "SaveFullCellSnapshot: electrolyte potential is null.");

    MFEM_VERIFY(
        state.RxnA_gf,
        "SaveFullCellSnapshot: anode reaction is null.");

    MFEM_VERIFY(
        state.RxnC_gf,
        "SaveFullCellSnapshot: cathode reaction is null.");

    MFEM_VERIFY(
        state.RxnE_gf,
        "SaveFullCellSnapshot: electrolyte reaction is null.");

    const std::string suffix = MakeTimestepSuffix(t);


    // ========================================================================
    // SAVE STATIC GEOMETRY AT t = 0
    // ========================================================================

    if (t == 0)
    {
        geometry.parallelMesh->SaveAsOne(
            (outdir + "/pmesh").c_str());

        domain_parameters.psiA->SaveAsOne(
            (outdir + "/psiA").c_str());

        domain_parameters.psiC->SaveAsOne(
            (outdir + "/psiC").c_str());

        domain_parameters.pse->SaveAsOne(
            (outdir + "/pse").c_str());


        // --------------------------------------------------------------------
        // Individual anode particle/material masks
        // --------------------------------------------------------------------

        for (int k = 0;
             k < static_cast<int>(domain_parameters.psA.size());
             ++k)
        {
            MFEM_VERIFY(
                domain_parameters.psA[k],
                "SaveFullCellSnapshot: anode particle phase field is null.");

            std::ostringstream name;

            name << outdir
                 << "/psA_"
                 << k;

            domain_parameters.psA[k]->SaveAsOne(
                name.str().c_str());
        }


        // --------------------------------------------------------------------
        // Individual cathode particle/material masks
        // --------------------------------------------------------------------

        for (int k = 0;
             k < static_cast<int>(domain_parameters.psC.size());
             ++k)
        {
            MFEM_VERIFY(
                domain_parameters.psC[k],
                "SaveFullCellSnapshot: cathode particle phase field is null.");

            std::ostringstream name;

            name << outdir
                 << "/psC_"
                 << k;

            domain_parameters.psC[k]->SaveAsOne(
                name.str().c_str());
        }
    }


    // ========================================================================
    // BUILD COMBINED ANODE AND CATHODE CONCENTRATIONS
    // ========================================================================

    mfem::ParGridFunction CnA(
        geometry.parfespace.get());

    mfem::ParGridFunction CnC(
        geometry.parfespace.get());


    BuildCombinedParticleConcentration(
        state.anode_particles,
        domain_parameters.psA,
        CnA);

    BuildCombinedParticleConcentration(
        state.cathode_particles,
        domain_parameters.psC,
        CnC);


    // ========================================================================
    // SAVE CONCENTRATIONS
    // ========================================================================

    CnA.SaveAsOne(
        (outdir + "/CnA" + suffix).c_str());

    CnC.SaveAsOne(
        (outdir + "/CnC" + suffix).c_str());

    state.CnE_gf->SaveAsOne(
        (outdir + "/CnE" + suffix).c_str());


    // ========================================================================
    // SAVE POTENTIALS
    // ========================================================================

    state.phA_gf->SaveAsOne(
        (outdir + "/phA" + suffix).c_str());

    state.phC_gf->SaveAsOne(
        (outdir + "/phC" + suffix).c_str());

    state.phE_gf->SaveAsOne(
        (outdir + "/phE" + suffix).c_str());


    // ========================================================================
    // SAVE REACTIONS
    // ========================================================================

    state.RxnA_gf->SaveAsOne(
        (outdir + "/RxnA" + suffix).c_str());

    state.RxnC_gf->SaveAsOne(
        (outdir + "/RxnC" + suffix).c_str());

    state.RxnE_gf->SaveAsOne(
        (outdir + "/RxnE" + suffix).c_str());
}