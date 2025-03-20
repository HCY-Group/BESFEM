/**
 * @file Reaction.cpp
 * @brief Implementation of the Reaction class for electrochemical reaction modeling in batteries.
 */

#include "Reaction.hpp"
#include "../code/Constants.hpp"
#include "mfem.hpp"
#include "CnE.hpp"
#include "PotE.hpp"
#include "PotP.hpp"


Reaction::Reaction(Initialize_Geometry &geo, Domain_Parameters &para)
    : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para), AvP(*para.AvP), AvB(*para.AvB), EVol(para.EVol)

{
    
    nE = geometry.nE; 
    nC = geometry.nC; 
    nV = geometry.nV; 

    // i0C = new mfem::ParGridFunction(fespace.get()); // exchange current density
    // OCV = new mfem::ParGridFunction(fespace.get()); // open circuit voltage
    // Kfw = new mfem::ParGridFunction(fespace.get()); // forward reaction constant
    // Kbw = new mfem::ParGridFunction(fespace.get()); // backward rection constant

    // dPHE = new mfem::ParGridFunction(fespace.get()); // voltage drop

    i0C = std::make_unique<mfem::ParGridFunction>(fespace.get());
    OCV = std::make_unique<mfem::ParGridFunction>(fespace.get());
    Kfw = std::make_unique<mfem::ParGridFunction>(fespace.get());
    Kbw = std::make_unique<mfem::ParGridFunction>(fespace.get());
    dPHE = std::make_unique<mfem::ParGridFunction>(fespace.get());

}

// Reaction::~Reaction()
// {
//     delete i0C;
//     delete OCV;
//     delete Kfw;
//     delete Kbw;
//     delete dPHE;
// }

void Reaction::Initialize(mfem::ParGridFunction &Rx, double initial_value) {

    SetInitialReaction(Rx, initial_value);
    Rx = AvP; // Scale by active particle surface area
    Rx *= 1.0e-8; // Apply a scaling factor

}

void Reaction::SetInitialReaction(mfem::ParGridFunction &Rx, double initial_value) {
    
    for (int i = 0; i < Rx.Size(); ++i) {
        Rx(i) = initial_value;
    }

}

void Reaction::ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2){

    for (int vi = 0; vi < nV; vi++){
        if ( AvB(vi) * Constants::dh > 0.0 ){ // Check for interface presence

            (*dPHE)(vi) = phx1(vi) - phx2(vi); // Voltage drop across the interface
            Rx(vi) = AvP(vi) * ((*Kfw)(vi)*Cn2(vi)*exp(-Constants::alp*Constants::Cst1*(*dPHE)(vi)) - \
					                   (*Kbw)(vi)*Cn1(vi)*exp( Constants::alp*Constants::Cst1*(*dPHE)(vi)));

        }
    }

}

// rate constants and exchange current density at interface
void Reaction::ExchangeCurrentDensity(mfem::ParGridFunction &Cn){

    for (int vi = 0; vi < nV; vi++){

        if(AvB(vi) * Constants::dh > 0.0){ 
            double val = -0.2 * (Cn(vi) - 0.37) - 1.559 - 0.9376 * tanh(8.961 * Cn(vi) - 3.195);
            (*i0C)(vi) = pow(10.0, val) * 1.0e-3; // Exchange current density

            (*OCV)(vi) = 1.095 * Cn(vi) * Cn(vi) - 8.324e-7 * exp(14.31 * Cn(vi)) + 4.692 * exp(-0.5389 * Cn(vi));

            (*Kfw)(vi) = (*i0C)(vi) / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * (*OCV)(vi));
            (*Kbw)(vi) = (*i0C)(vi) / (Constants::Frd * Cn(vi)) * exp(-Constants::alp * Constants::Cst1 * (*OCV)(vi));

        }
    }
}


void Reaction::TotalReactionCurrent(mfem::ParGridFunction &Rx, double &global_current){

    local_current = 0.0;
    mfem::Array<double> VtxVal(nC);
    mfem::Vector EAvg(nE);

    // Loop over all elements to calculate local reaction current
    for (int ei = 0; ei < nE; ei++){
        Rx.GetNodalValues(ei,VtxVal) ;
        // double val = 0.0;
        // for (int vt = 0; vt < nC; vt++){
        //     val += VtxVal[vt];
        // }
        double val = std::accumulate(VtxVal.begin(), VtxVal.end(), 0.0);
        EAvg(ei) = val/nC;	// Average reaction rate for the element
        local_current += EAvg(ei)*EVol(ei) ; // Accumulate local current		 			
    }
    
    // Perform global reduction to compute the total current
    MPI_Allreduce(&local_current, &global_current, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}