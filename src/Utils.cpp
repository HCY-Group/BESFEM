#include "../include/Utils.hpp"
#include <numeric>

Utils::Utils(Initialize_Geometry &geo, Domain_Parameters &para)
    : geometry_(geo), domain_(para), pmesh_(geo.parallelMesh.get()), fes_(geo.parfespace),
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

void Utils::SaveSimulationSnapshot(int t, const std::string &outdir, Initialize_Geometry &geometry, Domain_Parameters &domain_parameters, mfem::ParGridFunction &phA,
    mfem::ParGridFunction &phC, mfem::ParGridFunction &phE, mfem::ParGridFunction &CnA, mfem::ParGridFunction &CnC, mfem::ParGridFunction &CnE, mfem::ParGridFunction &CnApsi,
    mfem::ParGridFunction &CnCpsi, mfem::ParGridFunction &CnEpsi, mfem::ParGridFunction &CnP, int save_interval)
{
    if (t % save_interval != 0) return;

    std::ostringstream step;
    step << "_" << std::setw(5) << std::setfill('0') << t;
    std::string suff = step.str();

    geometry.parallelMesh->SaveAsOne((outdir + "/pmesh" + suff).c_str());

    domain_parameters.psi->SaveAsOne((outdir + "/psi" + suff).c_str());
    domain_parameters.pse->SaveAsOne((outdir + "/pse" + suff).c_str());

    phA.SaveAsOne((outdir + "/phA" + suff).c_str());
    phC.SaveAsOne((outdir + "/phC" + suff).c_str());
    phE.SaveAsOne((outdir + "/phE" + suff).c_str());

    CnA.SaveAsOne((outdir + "/CnA_raw" + suff).c_str());
    CnC.SaveAsOne((outdir + "/CnC_raw" + suff).c_str());
    CnE.SaveAsOne((outdir + "/CnE_raw" + suff).c_str());

    CnApsi = CnA; CnApsi *= *domain_parameters.psA;
    CnApsi.SaveAsOne((outdir + "/CnA" + suff).c_str());

    CnCpsi = CnC; CnCpsi *= *domain_parameters.psC;
    CnCpsi.SaveAsOne((outdir + "/CnC" + suff).c_str());

    CnEpsi = CnE; CnEpsi *= *domain_parameters.pse;
    CnEpsi.SaveAsOne((outdir + "/CnE" + suff).c_str());

    CnP = CnApsi;
    CnP += CnCpsi;
    CnP.SaveAsOne((outdir + "/CnP" + suff).c_str());
}

void Utils::SaveSimulationSnapshot(int t, const std::string &outdir,
    Initialize_Geometry &geometry, Domain_Parameters &domain_parameters, mfem::ParGridFunction &phC, mfem::ParGridFunction &phE, mfem::ParGridFunction &CnC, mfem::ParGridFunction &CnE,
    mfem::ParGridFunction &CnCpsi, mfem::ParGridFunction &CnEpsi, int save_interval)
{
    if (t % save_interval != 0) return;

    std::ostringstream step;
    step << "_" << std::setw(5) << std::setfill('0') << t;
    std::string suff = step.str();

    geometry.parallelMesh->SaveAsOne((outdir + "/pmesh" + suff).c_str());
    domain_parameters.psi->SaveAsOne((outdir + "/psi" + suff).c_str());
    domain_parameters.pse->SaveAsOne((outdir + "/pse" + suff).c_str());

    phC.SaveAsOne((outdir + "/phC" + suff).c_str());
    phE.SaveAsOne((outdir + "/phE" + suff).c_str());

    CnC.SaveAsOne((outdir + "/CnC_raw" + suff).c_str());
    CnE.SaveAsOne((outdir + "/CnE_raw" + suff).c_str());

    CnCpsi = CnC; CnCpsi *= *domain_parameters.psi;
    CnCpsi.SaveAsOne((outdir + "/CnC" + suff).c_str());

    CnEpsi = CnE; CnEpsi *= *domain_parameters.pse;
    CnEpsi.SaveAsOne((outdir + "/CnE" + suff).c_str());
}

void Utils::SetInitialValue(mfem::ParGridFunction &Cn, double initial_value)
    {
        for (int i = 0; i < Cn.Size(); i++)
            Cn(i) = initial_value;
    }

