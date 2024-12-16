/**
 * @file Reaction.cpp
 * @brief Implementation of the Reaction class for electrochemical reaction modeling in batteries.
 */

#include "Reaction.hpp"
#include "Mesh_Handler.hpp"
#include "mfem.hpp"
#include "CnE.hpp"
#include "PotE.hpp"
#include "PotP.hpp"


Reaction::Reaction(mfem::ParMesh *pm, mfem::ParFiniteElementSpace *fe, MeshHandler &mh)
    : pmesh(pm), fespace(fe), mesh_handler(mh), AvP(*mh.GetAvP()), AvB(*mh.GetAvB()), EVol(mh.GetElementVolume())

{
    
    nE = mesh_handler.GetNE(); // Number of elements.
    nC = mesh_handler.GetNC(); // Number of corners.
    nV = mesh_handler.GetNV(); // Number of vertices.

    i0C = new mfem::ParGridFunction(fespace); // exchange current density
    OCV = new mfem::ParGridFunction(fespace); // open circuit voltage
    Kfw = new mfem::ParGridFunction(fespace); // forward reaction constant
    Kbw = new mfem::ParGridFunction(fespace); // backward rection constant

    dPHE = new mfem::ParGridFunction(fespace); // voltage drop



}

void Reaction::Initialize(mfem::ParGridFunction &Rx, double initial_value) {

    SetInitialReaction(Rx, initial_value);
    Rx = AvP; // Scale by active particle surface area.
    Rx *= 1.0e-8; // Apply a scaling factor.

}

void Reaction::SetInitialReaction(mfem::ParGridFunction &Rx, double initial_value) {
    
    for (int i = 0; i < Rx.Size(); ++i) {
        Rx(i) = initial_value;
    }

}

void Reaction::ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2){

    for (int vi = 0; vi < nV; vi++){
        if ( AvB(vi) * Constants::dh > 0.0 ){ // Check for interface presence.

            (*dPHE)(vi) = phx1(vi) - phx2(vi); // Voltage drop across the interface.
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
            (*i0C)(vi) = pow(10.0, val) * 1.0e-3; // Exchange current density.

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

    // Loop over all elements to calculate local reaction current.
    for (int ei = 0; ei < nE; ei++){
        Rx.GetNodalValues(ei,VtxVal) ;
        double val = 0.0;
        for (int vt = 0; vt < nC; vt++){
            val += VtxVal[vt];
        }
        EAvg(ei) = val/nC;	// Average reaction rate for the element.
        local_current += EAvg(ei)*EVol(ei) ; // Accumulate local current.		 			
    }
    
    // Perform global reduction to compute the total current.
    MPI_Allreduce(&local_current, &global_current, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}