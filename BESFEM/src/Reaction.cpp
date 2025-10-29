/**
 * @file Reaction.cpp
 * @brief Implementation of the Reaction class for electrochemical reaction modeling in batteries.
 */

 #include "../include/Reaction.hpp"
 #include "../inputs/Constants.hpp"
 #include "mfem.hpp"
 #include "../include/CnE.hpp"
 #include "../include/PotE.hpp"
 #include "../include/PotC.hpp"
#include "../include/PotA.hpp"
 #include <fstream>
 #include <cmath>
 
 Reaction::Reaction(Initialize_Geometry &geo, Domain_Parameters &para)
     : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo),
       domain_parameters(para), EVol(para.EVol),
    //    AvP(para.AvP ? para.AvP.get() : nullptr), AvB(para.AvB ? para.AvB.get() : nullptr),
    //    AvA(para.AvA ? para.AvA.get() : nullptr), AvC(para.AvC ? para.AvC.get() : nullptr)
    AvP(para.AvP ? para.AvP.get() : nullptr),
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

 
    std::ifstream myXfile("../inputs/C_Li_X_101.txt"); // ticks
    std::ifstream mydFfile("../inputs/C_Li_M6_101.txt"); // chemical potential
    std::ifstream myMBfile("../inputs/C_Li_Mb5_101.txt"); // mobility
    std::ifstream myOCVfile("../inputs/C_Li_O3_101.txt"); // OCV
    std::ifstream myi0file("../inputs/C_Li_J2_101.txt"); // i0
    
    for (int i = 0; i < 101; i++) myXfile >> Ticks(i);
    for (int i = 0; i < 101; i++) mydFfile >> chmPot(i);
    for (int i = 0; i < 101; i++) myMBfile >> Mobility(i);
    for (int i = 0; i < 101; i++) myOCVfile >> OCV_file(i);
    for (int i = 0; i < 101; i++) myi0file >> i0_file(i);
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

// rate constants and exchange current density at interface
void Reaction::ExchangeCurrentDensity(mfem::ParGridFunction &Cn){
    for (int vi = 0; vi < nV; vi++){
        if((*AvB)(vi) * Constants::dh > 0.0){ 
            double val = -0.2 * (Cn(vi) - 0.37) - 1.559 - 0.9376 * tanh(8.961 * Cn(vi) - 3.195);
            (*i0C)(vi) = pow(10.0, val) * 1.0e-3; // Exchange current density
            (*OCV)(vi) = 1.095 * Cn(vi) * Cn(vi) - 8.324e-7 * exp(14.31 * Cn(vi)) + 4.692 * exp(-0.5389 * Cn(vi)); // open circuit voltage
            (*Kfw)(vi) = (*i0C)(vi) / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * (*OCV)(vi)); // forward reaction constant
            (*Kbw)(vi) = (*i0C)(vi) / (Constants::Frd * Cn(vi)) * exp(-Constants::alp * Constants::Cst1 * (*OCV)(vi)); // backward rection constant
        }
    }
}

void Reaction::ExchangeCurrentDensity(mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2){
    for (int vi = 0; vi < nV; vi++){
        if((*AvC)(vi) * Constants::dh > 0.0){ 
            double val = -0.2 * (Cn1(vi) - 0.37) - 1.559 - 0.9376 * tanh(8.961 * Cn1(vi) - 3.195);
            (*i0CC)(vi) = pow(10.0, val) * 1.0e-3; // Exchange current density
            (*OCVC)(vi) = 1.095 * Cn1(vi) * Cn1(vi) - 8.324e-7 * exp(14.31 * Cn1(vi)) + 4.692 * exp(-0.5389 * Cn1(vi)); // open circuit voltage
            (*KfC)(vi) = (*i0CC)(vi) / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * (*OCVC)(vi)); // forward reaction constant
            (*KbC)(vi) = (*i0CC)(vi) / (Constants::Frd * Cn1(vi)) * exp(-Constants::alp * Constants::Cst1 * (*OCVC)(vi)); // backward rection constant
        }

        if((*AvA)(vi) * Constants::dh > 0.0){
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

void Reaction::TableExchangeCurrentDensity(mfem::ParGridFunction &Cn)
{
    for (int vi = 0; vi < nV; vi++) {
        double cn_val = Cn(vi);

        double i0 = GetTableValues(cn_val, Ticks, i0_file) * 1.0e-3; // Convert mA to A
        double ocv = GetTableValues(cn_val, Ticks, OCV_file);

        (*i0C)(vi) = i0;
        (*OCV)(vi) = ocv;
        (*Kfw)(vi) = i0 / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * ocv);
        (*Kbw)(vi) = i0 / (Constants::Frd * cn_val) * exp(-Constants::alp * Constants::Cst1 * ocv);
    
    }
}


void Reaction::ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2)
{
    for (int vi = 0; vi < nV; vi++){
        if ( (*AvB)(vi) * Constants::dh > 0.0 ){ // Check for interface presence
            (*dPHE)(vi) = phx1(vi) - phx2(vi); // Voltage drop across the interface
            Rx(vi) = (*AvP)(vi) * ((*Kfw)(vi)*Cn2(vi)*exp(-Constants::alp*Constants::Cst1*(*dPHE)(vi)) - \
					                   (*Kbw)(vi)*Cn1(vi)*exp( Constants::alp*Constants::Cst1*(*dPHE)(vi)));

        }
    }
}

void Reaction::ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &Cn3, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2, mfem::ParGridFunction &phx3)
{
    for (int vi = 0; vi < nV; vi++){
        
        if ( (*AvA)(vi) * Constants::dh > 0.0 ){ // Check for interface presence
                (*dPHA)(vi) = phx2(vi) - phx3(vi); // Voltage drop across the interface
                Rx2(vi) = (*AvA)(vi) * ((*KfA)(vi)*Cn3(vi)*exp(-Constants::alp*Constants::Cst1*(*dPHA)(vi)) - \
                                           (*KbA)(vi)*Cn2(vi)*exp( Constants::alp*Constants::Cst1*(*dPHA)(vi)));
            }

        if ( (*AvC)(vi) * Constants::dh > 0.0 ){ // Check for interface presence
                (*dPHC)(vi) = phx1(vi) - phx3(vi); // Voltage drop across the interface
                Rx1(vi) = (*AvC)(vi) * ((*KfC)(vi)*Cn3(vi)*exp(-Constants::alp*Constants::Cst1*(*dPHC)(vi)) - \
                                           (*KbC)(vi)*Cn1(vi)*exp( Constants::alp*Constants::Cst1*(*dPHC)(vi)));
            }
    }

    // std::cout << "dPHC min/max: " << dPHC->Min() << " " << dPHC->Max() << std::endl;
    // std::cout << "dPHA min/max: " << dPHA->Min() << " " << dPHA->Max() << std::endl;
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

    //  if (mfem::Mpi::WorldRank() == 0){std::cout << "local current: " << local_current << std::endl;}
 
     MPI_Allreduce(&local_current, &global_current, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 }