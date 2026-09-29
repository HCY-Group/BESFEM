#include "../include/Reaction.hpp"
#include "../include/Constants.hpp"
#include "mfem.hpp"
#include "../include/MaterialProperties.hpp"
#include "../include/SimulationConfig.hpp"
#include "../include/SimTypes.hpp"
#include <fstream>
#include <cmath>
 
Reaction::Reaction(Initialize_Geometry &geo, Domain_Parameters &para, const SimulationConfig &cfg)
    : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), cfg(cfg),
    domain_parameters(para), EVol(para.EVol)
{
    nE = geometry.nE; 
    nC = geometry.nC; 
    nV = geometry.nV; 

    i0C = std::make_unique<mfem::ParGridFunction>(fespace.get()); // exchange current density
    OCV = std::make_unique<mfem::ParGridFunction>(fespace.get()); // open circuit voltage

    Kfw = std::make_unique<mfem::ParGridFunction>(fespace.get()); // forward reaction constant
    Kbw = std::make_unique<mfem::ParGridFunction>(fespace.get()); // backward reaction constant

    dPHE = std::make_unique<mfem::ParGridFunction>(fespace.get()); // voltage drop

}

double Reaction::GetTableValues(double cn, const mfem::Vector &ticks, const mfem::Vector &data)
{
    if (cn < 1.0e-6) cn = 1.0e-6;
    if (cn > 0.999999) cn = 0.999999;

    int idx = std::floor(cn / 0.01);
    if (idx < 0) idx = 0;
    if (idx > 99) idx = 99;

    return data(idx) + (cn - ticks(idx)) / 0.01 * (data(idx + 1) - data(idx));
}

void Reaction::Initialize(mfem::ParGridFunction &Rx, double initial_value)
{
    SetInitialReaction(Rx, initial_value);
}

void Reaction::SetInitialReaction(mfem::ParGridFunction &Rx, double initial_value)
{
    for (int i = 0; i < Rx.Size(); ++i) {
        Rx(i) = initial_value;
}
}

void Reaction::ExchangeCurrentDensity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &AvP_in, sim::MaterialType material)
{
    static bool printed = false;

    *i0C = 0.0;
    *OCV = 0.0;
    *Kfw = 0.0;
    *Kbw = 0.0;

    for (int vi = 0; vi < nV; vi++)
    {
        if ((AvP_in)(vi) * cfg.dh > 1e-3)
        {
            const double cn_val = Cn(vi);

            (*i0C)(vi) = MaterialProperties::ExchangeCurrentDensity(material, cn_val);
            (*OCV)(vi) = MaterialProperties::OCV(material, cn_val);

            (*Kfw)(vi) = (*i0C)(vi) / (Constants::Frd * 0.001) * std::exp(Constants::alp * Constants::Cst1 * (*OCV)(vi));
            (*Kbw)(vi) = (*i0C)(vi) / (Constants::Frd * cn_val) * std::exp(-Constants::alp * Constants::Cst1 * (*OCV)(vi));
        }
    }
}

void Reaction::TableExchangeCurrentDensity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &AvP_in)
{
    for (int vi = 0; vi < nV; vi++) {

        if ((AvP_in)(vi) * cfg.dh > 1e-3) { // Check for interface presence
            double cn_val = Cn(vi);

            double i0 = GetTableValues(cn_val, Ticks, i0_file) * 1.0e-3; // Convert mA to A
            double ocv = GetTableValues(cn_val, Ticks, OCV_file);

            (*i0C)(vi) = i0;
            (*OCV)(vi) = ocv;
            (*Kfw)(vi) = i0 / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * ocv);
            (*Kbw)(vi) = i0 / (Constants::Frd * cn_val) * exp(-Constants::alp * Constants::Cst1 * ocv);
        }

    }
}

void Reaction::ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2, mfem::ParGridFunction &AvP_in)
{
    Rx = 0.0;

    for (int vi = 0; vi < nV; vi++){
        if ( (AvP_in)(vi) * cfg.dh > 1e-3){ // Check for interface presence
            (*dPHE)(vi) = phx1(vi) - phx2(vi); // Voltage drop across the interface
            Rx(vi) = (AvP_in)(vi) * ((*Kfw)(vi)*Cn2(vi)*exp(-Constants::alp*Constants::Cst1*(*dPHE)(vi)) - \
                                        (*Kbw)(vi)*Cn1(vi)*exp( Constants::alp*Constants::Cst1*(*dPHE)(vi)));

        }
    }
}

void Reaction::TotalReactionCurrent(mfem::ParGridFunction &Rx, double &global_current)
{
    local_current = 0.0;
    mfem::Array<double> VtxVal(nC);
    mfem::Vector EAvg(nE);

    for (int ei = 0; ei < nE; ++ei) {
        Rx.GetNodalValues(ei, VtxVal);
        double sum = std::accumulate(VtxVal.begin(), VtxVal.end(), 0.0);
        EAvg(ei) = sum / nC;
        local_current += EAvg(ei) * EVol(ei);
    }

    MPI_Allreduce(&local_current, &global_current, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
