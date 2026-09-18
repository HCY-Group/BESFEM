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
    domain_parameters(para), EVol(para.EVol),
AvB(para.AvB ? para.AvB.get() : nullptr),
AvA(para.AvA ? para.AvA.get() : nullptr),
AvC(para.AvC ? para.AvC.get() : nullptr)
{
nE = geometry.nE; 
nC = geometry.nC; 
nV = geometry.nV; 

i0C = std::make_unique<mfem::ParGridFunction>(fespace.get()); // exchange current density
OCV = std::make_unique<mfem::ParGridFunction>(fespace.get()); // open circuit voltage

i0CC = std::make_unique<mfem::ParGridFunction>(fespace.get()); // exchange current density (cathode)
OCVC = std::make_unique<mfem::ParGridFunction>(fespace.get()); // open circuit voltage (cathode)
i0CA = std::make_unique<mfem::ParGridFunction>(fespace.get()); // exchange current density (anode)
OCVA = std::make_unique<mfem::ParGridFunction>(fespace.get()); // open circuit voltage (anode)

Kfw = std::make_unique<mfem::ParGridFunction>(fespace.get()); // forward reaction constant
Kbw = std::make_unique<mfem::ParGridFunction>(fespace.get()); // backward reaction constant

KfA = std::make_unique<mfem::ParGridFunction>(fespace.get()); // forward reaction constant (anode)
KbA = std::make_unique<mfem::ParGridFunction>(fespace.get()); // backward reaction constant (anode)

KfC = std::make_unique<mfem::ParGridFunction>(fespace.get()); // forward reaction constant (cathode)
KbC = std::make_unique<mfem::ParGridFunction>(fespace.get()); // backward reaction constant (cathode)

dPHE = std::make_unique<mfem::ParGridFunction>(fespace.get()); // voltage drop
dPHA = std::make_unique<mfem::ParGridFunction>(fespace.get()); // voltage drop
dPHC = std::make_unique<mfem::ParGridFunction>(fespace.get()); // voltage drop

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

void Reaction::ExchangeCurrentDensity(mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &AvA_in, mfem::ParGridFunction &AvC_in){
    for (int vi = 0; vi < nV; vi++){
        if((AvC_in)(vi) * cfg.dh > 1e-3){ 
            double val = -0.2 * (Cn1(vi) - 0.37) - 1.559 - 0.9376 * tanh(8.961 * Cn1(vi) - 3.195);
            (*i0CC)(vi) = pow(10.0, val) * 1.0e-3; // Exchange current density
            (*OCVC)(vi) = 1.095 * Cn1(vi) * Cn1(vi) - 8.234e-7 * exp(14.31 * Cn1(vi)) + 4.692 * exp(-0.5389 * Cn1(vi)); // open circuit voltage
            (*KfC)(vi) = (*i0CC)(vi) / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * (*OCVC)(vi)); // forward reaction constant
            (*KbC)(vi) = (*i0CC)(vi) / (Constants::Frd * Cn1(vi)) * exp(-Constants::alp * Constants::Cst1 * (*OCVC)(vi)); // backward rection constant
        }

        if((AvA_in)(vi) * cfg.dh > 1e-3){
            double cn_val = Cn2(vi);
            double i0 = GetTableValues(cn_val, Ticks, i0_file) * 1.0e-3; // Convert mA to A
            double ocv = GetTableValues(cn_val, Ticks, OCV_file);

            (*i0CA)(vi) = i0;
            (*OCVA)(vi) = ocv;
            (*KfA)(vi) = i0 / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * ocv);
            (*KbA)(vi) = i0 / (Constants::Frd * cn_val) * exp(-Constants::alp * Constants::Cst1 * ocv);
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

void Reaction::ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &Cn3, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2, mfem::ParGridFunction &phx3, mfem::ParGridFunction &AvA_in, mfem::ParGridFunction &AvC_in)
{
    for (int vi = 0; vi < nV; vi++){
        
        if ( (AvA_in)(vi) * cfg.dh > 1e-3 ){ // Check for interface presence
                (*dPHA)(vi) = phx2(vi) - phx3(vi); // Voltage drop across the interface
                Rx2(vi) = (AvA_in)(vi) * ((*KfA)(vi)*Cn3(vi)*exp(-Constants::alp*Constants::Cst1*(*dPHA)(vi)) - \
                                            (*KbA)(vi)*Cn2(vi)*exp( Constants::alp*Constants::Cst1*(*dPHA)(vi)));
            }

        if ( (AvC_in)(vi) * cfg.dh > 1e-3 ){ // Check for interface presence
                (*dPHC)(vi) = phx1(vi) - phx3(vi); // Voltage drop across the interface
                Rx1(vi) = (AvC_in)(vi) * ((*KfC)(vi)*Cn3(vi)*exp(-Constants::alp*Constants::Cst1*(*dPHC)(vi)) - \
                                            (*KbC)(vi)*Cn1(vi)*exp( Constants::alp*Constants::Cst1*(*dPHC)(vi)));
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