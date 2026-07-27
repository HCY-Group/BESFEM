#include "../include/Utils.hpp"
#include "../include/Constants.hpp"
#include "../include/SimulationConfig.hpp"
#include "../include/MaterialProperties.hpp"

#include <numeric>

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

void Utils::CalculateReactionInfx(mfem::ParGridFunction &Rx, double &xCrnt)
{
    xCrnt = 0.0;

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

void Utils::ComputePairFlux(mfem::ParGridFunction &sum_part, mfem::ParGridFunction &weight, mfem::ParGridFunction &grad_psi, mfem::ParGridFunction &mu_1, mfem::ParGridFunction &mu_2)
{
    for (int vi = 0; vi < nV_; vi++){

        double grad_psi_val = grad_psi(vi);
        double weight_val = weight(vi);
        double mu1_val = mu_1(vi);
        double mu2_val = mu_2(vi);

        const double rho = MaterialProperties::SiteDensity(cfg.cathode_materials[0]);

        sum_part(vi) = weight_val * grad_psi_val * rho * (1.0/Constants::RT) * Constants::Perm * (mu2_val - mu1_val);
    }

}

// Full Cell
void Utils::SaveSimulationSnapshot(int t, const std::string &outdir, Initialize_Geometry &geometry, Domain_Parameters &domain_parameters, mfem::ParGridFunction &phA,
    mfem::ParGridFunction &phC, mfem::ParGridFunction &phE, mfem::ParGridFunction &CnA, mfem::ParGridFunction &CnC, mfem::ParGridFunction &CnE, mfem::ParGridFunction &CnApsi,
    mfem::ParGridFunction &CnCpsi, mfem::ParGridFunction &CnEpsi, mfem::ParGridFunction &CnP, int save_interval)
{
    if (t % save_interval != 0) return;

    std::ostringstream step;
    step << "_" << std::setw(5) << std::setfill('0') << t;
    std::string suff = step.str();

    if (t == 0)
    {
        geometry.parallelMesh->SaveAsOne((outdir + "/pmesh").c_str());
        domain_parameters.psiA->SaveAsOne((outdir + "/psiA").c_str());
        domain_parameters.psiC->SaveAsOne((outdir + "/psiC").c_str());
        domain_parameters.pse->SaveAsOne((outdir + "/pse").c_str());
    }

    phA.SaveAsOne((outdir + "/phA" + suff).c_str());
    phC.SaveAsOne((outdir + "/phC" + suff).c_str());
    phE.SaveAsOne((outdir + "/phE" + suff).c_str());

    CnA.SaveAsOne((outdir + "/CnA_raw" + suff).c_str());
    CnC.SaveAsOne((outdir + "/CnC_raw" + suff).c_str());
    CnE.SaveAsOne((outdir + "/CnE_raw" + suff).c_str());

    CnApsi = CnA; CnApsi *= *domain_parameters.psiA;
    CnApsi.SaveAsOne((outdir + "/CnA" + suff).c_str());

    CnCpsi = CnC; CnCpsi *= *domain_parameters.psiC;
    CnCpsi.SaveAsOne((outdir + "/CnC" + suff).c_str());

    CnEpsi = CnE; CnEpsi *= *domain_parameters.pse;
    CnEpsi.SaveAsOne((outdir + "/CnE" + suff).c_str());

    CnP = CnApsi;
    CnP += CnCpsi;
    CnP.SaveAsOne((outdir + "/CnP" + suff).c_str());
}

// Half Cell
void Utils::SaveSimulationSnapshot(int t, const std::string &outdir,
    Initialize_Geometry &geometry, Domain_Parameters &domain_parameters, mfem::ParGridFunction &phC, mfem::ParGridFunction &phE, mfem::ParGridFunction &CnC, mfem::ParGridFunction &CnE,
    mfem::ParGridFunction &CnCpsi, mfem::ParGridFunction &CnEpsi, int save_interval)
{
    if (t % save_interval != 0) return;

    std::ostringstream step;
    step << "_" << std::setw(5) << std::setfill('0') << t;
    std::string suff = step.str();

    if (t == 0)
    {
        geometry.parallelMesh->SaveAsOne((outdir + "/pmesh").c_str());
        domain_parameters.psi->SaveAsOne((outdir + "/psi").c_str());
        domain_parameters.pse->SaveAsOne((outdir + "/pse").c_str());
    }

    phC.SaveAsOne((outdir + "/phC" + suff).c_str());
    phE.SaveAsOne((outdir + "/phE" + suff).c_str());

    CnC.SaveAsOne((outdir + "/CnP_raw" + suff).c_str());
    CnE.SaveAsOne((outdir + "/CnE_raw" + suff).c_str());

    CnCpsi = CnC; CnCpsi *= *domain_parameters.psi;
    CnCpsi.SaveAsOne((outdir + "/CnP" + suff).c_str());

    CnEpsi = CnE; CnEpsi *= *domain_parameters.pse;
    CnEpsi.SaveAsOne((outdir + "/CnE" + suff).c_str());
}

void Utils::SetInitialValue(mfem::ParGridFunction &Cn, double initial_value)
    {
        for (int i = 0; i < Cn.Size(); i++)
            Cn(i) = initial_value;
    }

// void Utils::SaveSimulationSnapshotMulti(int t, const std::string &outdir, Initialize_Geometry &geometry,
//     Domain_Parameters &domain_parameters, const std::vector<mfem::ParGridFunction*> &particle_cn,
//     std::vector<std::unique_ptr<mfem::ParGridFunction>> &particle_out, int save_interval)

//     {
//         if (t % save_interval != 0) return;

//         const int np = static_cast<int>(particle_cn.size());
        
//         std::ostringstream step;
//         step << "_" << std::setw(5) << std::setfill('0') << t;
//         const std::string suff = step.str();

//         if (t == 0)
//         {
//             geometry.parallelMesh->SaveAsOne((outdir + "/pmesh").c_str());
//             domain_parameters.psi->SaveAsOne((outdir + "/psi").c_str());
//             domain_parameters.pse->SaveAsOne((outdir + "/pse").c_str());
//         }

//         // Save each particle concentration and masked version
//         for (int k = 0; k < np; ++k)
//         {
//             std::ostringstream raw_name, masked_name;
//             raw_name << outdir << "/CnC_" << (k + 1) << suff;
//             masked_name << outdir << "/C" << (k + 1) << "_out" << suff;

//             // particle_cn[k]->SaveAsOne(raw_name.str().c_str());

//             *particle_out[k] = *particle_cn[k];
//             *particle_out[k] *= *domain_parameters.ps[k];
//             // particle_out[k]->SaveAsOne(masked_name.str().c_str());
//         }

//         // Build union mask and denominator
//         mfem::ParGridFunction psi_union(geometry.parfespace.get());
//         mfem::ParGridFunction denom(geometry.parfespace.get());
//         mfem::ParGridFunction CnP_total(geometry.parfespace.get());

//         psi_union = 0.0;
//         denom = 0.0;
//         CnP_total = 0.0;

//         for (int k = 0; k < np; ++k)
//         {
//             psi_union += *domain_parameters.ps[k];
//             denom += *domain_parameters.ps[k];
//         }

//         for (int i = 0; i < psi_union.Size(); ++i)
//         {
//             psi_union(i) = std::min(1.0, psi_union(i));
//         }

//         const double eps = 1e-30;

//         for (int i = 0; i < denom.Size(); ++i)
//         {
//             const double d = denom(i);

//             if (d > eps)
//             {
//                 double num = 0.0;
//                 for (int k = 0; k < np; ++k)
//                 {
//                     num += (*domain_parameters.ps[k])(i) * (*particle_cn[k])(i);
//                 }
//                 CnP_total(i) = num / (d + eps);
//             }
//             else
//             {
//                 CnP_total(i) = 0.0;
//             }
//         }

//         CnP_total *= psi_union;
//         CnP_total.SaveAsOne((outdir + "/CnP_total" + suff).c_str());
//     }

void Utils::SaveSimulationSnapshotMulti(int t, const std::string& outdir, Initialize_Geometry& geometry,
    Domain_Parameters& domain_parameters, const std::vector<mfem::ParGridFunction*>& particle_cn,
    const std::vector<std::unique_ptr<mfem::ParGridFunction>>& particle_ps, mfem::ParGridFunction& electrode_psi,
    std::vector<std::unique_ptr<mfem::ParGridFunction>>& particle_out, const std::string& electrode_name, int save_interval)
    {
        if (t % save_interval != 0)
        {
            return;
        }

        const int np = static_cast<int>(particle_cn.size());
        MFEM_VERIFY(static_cast<int>(particle_ps.size()) == np, "SaveSimulationSnapshotMulti: particle concentration and phase-field counts differ.");
        MFEM_VERIFY(static_cast<int>(particle_out.size()) == np, "SaveSimulationSnapshotMulti: particle concentration and output counts differ.");

        std::ostringstream step;
        step << "_" << std::setw(5) << std::setfill('0') << t;
        const std::string suff = step.str();

        // =====================================================================
        // Save mesh and phase fields at timestep zero
        // =====================================================================

        if (t == 0)
        {
            MFEM_VERIFY(geometry.parallelMesh, "SaveSimulationSnapshotMulti: parallel mesh is null.");
            MFEM_VERIFY(domain_parameters.pse, "SaveSimulationSnapshotMulti: electrolyte phase field is null.");

            geometry.parallelMesh->SaveAsOne( (outdir + "/pmesh").c_str());
            electrode_psi.SaveAsOne((outdir + "/psi" + electrode_name).c_str());
            domain_parameters.pse->SaveAsOne((outdir + "/pse").c_str());

            for (int k = 0; k < np; ++k)
            {
                MFEM_VERIFY(particle_ps[k], "SaveSimulationSnapshotMulti: null particle phase field.");
                std::ostringstream psi_name;
                psi_name << outdir << "/ps" << electrode_name << "_" << (k + 1);
                particle_ps[k]->SaveAsOne(psi_name.str().c_str());
            }
        }

        // =====================================================================
        // Save each particle concentration and masked concentration
        // =====================================================================

        for (int k = 0; k < np; ++k)
        {
            MFEM_VERIFY(particle_cn[k] != nullptr, "SaveSimulationSnapshotMulti: null particle concentration field.");
            MFEM_VERIFY(particle_ps[k], "SaveSimulationSnapshotMulti: null particle phase field.");
            MFEM_VERIFY(particle_out[k], "SaveSimulationSnapshotMulti: null particle output field.");
            MFEM_VERIFY(particle_cn[k]->Size() == particle_ps[k]->Size(), "SaveSimulationSnapshotMulti: concentration and phase-field sizes differ.");
            MFEM_VERIFY(particle_out[k]->Size() == particle_cn[k]->Size(), "SaveSimulationSnapshotMulti: output and concentration sizes differ.");

            std::ostringstream raw_name;
            std::ostringstream masked_name;

            // raw_name << outdir << "/Cn" << electrode_name << "_" << (k + 1) << suff;
            // masked_name << outdir << "/C" << electrode_name << "_" << (k + 1) << "_out" << suff;

            // particle_cn[k]->SaveAsOne(raw_name.str().c_str());
            *particle_out[k] = *particle_cn[k];
            *particle_out[k] *= *particle_ps[k];
            // particle_out[k]->SaveAsOne(masked_name.str().c_str());
        }

        // =====================================================================
        // Build union mask and total particle concentration
        // =====================================================================

        mfem::ParGridFunction psi_union(geometry.parfespace.get());
        mfem::ParGridFunction denom(geometry.parfespace.get());
        mfem::ParGridFunction CnP_total(geometry.parfespace.get());

        psi_union = 0.0;
        denom = 0.0;
        CnP_total = 0.0;

        for (int k = 0; k < np; ++k)
        {
            psi_union += *particle_ps[k];
            denom += *particle_ps[k];
        }

        for (int i = 0; i < psi_union.Size(); ++i)
        {
            psi_union(i) = std::min(1.0, psi_union(i));
        }

        const double eps = 1.0e-30;

        for (int i = 0; i < denom.Size(); ++i)
        {
            const double d = denom(i);

            if (d > eps)
            {
                double numerator = 0.0;

                for (int k = 0; k < np; ++k)
                {
                    numerator += (*particle_ps[k])(i) * (*particle_cn[k])(i);
                }

                CnP_total(i) = numerator / (d + eps);
            }
            else
            {
                CnP_total(i) = 0.0;
            }
        }

        CnP_total *= psi_union;
        CnP_total.SaveAsOne((outdir + "/Cn" + electrode_name + "_total" + suff).c_str());
    }

void Utils::SaveCombinedElectrodeSnapshot(
    int t,
    const std::string& outdir,
    Initialize_Geometry& geometry,
    const std::vector<mfem::ParGridFunction*>& anode_cn,
    const std::vector<std::unique_ptr<mfem::ParGridFunction>>& anode_ps,
    const std::vector<mfem::ParGridFunction*>& cathode_cn,
    const std::vector<std::unique_ptr<mfem::ParGridFunction>>& cathode_ps,
    int save_interval)
{
    if (t % save_interval != 0)
    {
        return;
    }

    MFEM_VERIFY(
        anode_cn.size() == anode_ps.size(),
        "SaveCombinedElectrodeSnapshot: anode concentration and phase-field counts differ.");

    MFEM_VERIFY(
        cathode_cn.size() == cathode_ps.size(),
        "SaveCombinedElectrodeSnapshot: cathode concentration and phase-field counts differ.");

    MFEM_VERIFY(
        geometry.parfespace,
        "SaveCombinedElectrodeSnapshot: finite-element space is null.");

    std::ostringstream step;
    step << "_" << std::setw(5) << std::setfill('0') << t;
    const std::string suffix = step.str();

    mfem::ParGridFunction electrode_union(geometry.parfespace.get());
    mfem::ParGridFunction numerator(geometry.parfespace.get());
    mfem::ParGridFunction denominator(geometry.parfespace.get());
    mfem::ParGridFunction combined_concentration(geometry.parfespace.get());

    electrode_union = 0.0;
    numerator = 0.0;
    denominator = 0.0;
    combined_concentration = 0.0;

    // =====================================================================
    // Add all anode particles
    // =====================================================================

    for (std::size_t k = 0; k < anode_cn.size(); ++k)
    {
        MFEM_VERIFY(
            anode_cn[k] != nullptr,
            "SaveCombinedElectrodeSnapshot: null anode concentration field.");

        MFEM_VERIFY(
            anode_ps[k] != nullptr,
            "SaveCombinedElectrodeSnapshot: null anode particle phase field.");

        MFEM_VERIFY(
            anode_cn[k]->Size() == anode_ps[k]->Size(),
            "SaveCombinedElectrodeSnapshot: anode concentration and phase-field sizes differ.");

        for (int i = 0; i < numerator.Size(); ++i)
        {
            const double psi_value = (*anode_ps[k])(i);

            numerator(i) += psi_value * (*anode_cn[k])(i);
            denominator(i) += psi_value;
        }
    }

    // =====================================================================
    // Add all cathode particles
    // =====================================================================

    for (std::size_t k = 0; k < cathode_cn.size(); ++k)
    {
        MFEM_VERIFY(
            cathode_cn[k] != nullptr,
            "SaveCombinedElectrodeSnapshot: null cathode concentration field.");

        MFEM_VERIFY(
            cathode_ps[k] != nullptr,
            "SaveCombinedElectrodeSnapshot: null cathode particle phase field.");

        MFEM_VERIFY(
            cathode_cn[k]->Size() == cathode_ps[k]->Size(),
            "SaveCombinedElectrodeSnapshot: cathode concentration and phase-field sizes differ.");

        for (int i = 0; i < numerator.Size(); ++i)
        {
            const double psi_value = (*cathode_ps[k])(i);

            numerator(i) += psi_value * (*cathode_cn[k])(i);
            denominator(i) += psi_value;
        }
    }

    // =====================================================================
    // Construct combined electrode mask and concentration
    // =====================================================================

    const double eps = 1.0e-30;

    for (int i = 0; i < combined_concentration.Size(); ++i)
    {
        electrode_union(i) = std::min(1.0, denominator(i));

        if (denominator(i) > eps)
        {
            combined_concentration(i) =
                numerator(i) / denominator(i);
        }
        else
        {
            combined_concentration(i) = 0.0;
        }
    }

    // Make the value exactly zero outside both electrodes.
    combined_concentration *= electrode_union;

    combined_concentration.SaveAsOne(
        (outdir + "/CnElectrodes_total" + suffix).c_str());

    // The mask only needs to be saved once because it does not vary in time.
    if (t == 0)
    {
        electrode_union.SaveAsOne(
            (outdir + "/psiElectrodes").c_str());
    }
}