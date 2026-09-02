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