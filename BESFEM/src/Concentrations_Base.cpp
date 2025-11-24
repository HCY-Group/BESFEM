

#include "../include/Concentrations_Base.hpp"
#include "../inputs/Constants.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "mfem.hpp"

ConcentrationBase::ConcentrationBase(Initialize_Geometry &geo, Domain_Parameters &para)
    : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para), EVol(para.EVol), gtPsi(para.gtPsi), gtPse(para.gtPse),
    nE(geo.nE), nC(geo.nC), nV(geo.nV), VtxVal(geo.nC), EAvg(geo.nE), gmesh(geo.globalMesh.get())


{
    // Allocate once for reuse
    VtxVal.SetSize(nC); // Set size for nodal values
    EAvg.SetSize(nE);   // Set size for average element contributions
    TmpF = std::make_unique<mfem::ParGridFunction>(fespace.get());

}

// void ConcentrationBase::SetupField(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);
// void ConcentrationBase::UpdateConcentration(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx);

// // void Concentrations::SetInitialConcentration(mfem::ParGridFunction &Cn, double initial_value) {
    
// //     for (int i = 0; i < Cn.Size(); ++i) {
// //         Cn(i) = initial_value;
// //     }

// // }

// // void Concentrations::LithiationCalculation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, double gtps) {
    
// //     // Temporary grid function to store the product of concentration and potential
// //     *TmpF = Cn; // Copy concentration values to the temporary grid function
// //     *TmpF *= psx; // Element-wise multiply concentration by potential

// //     double lSum = 0.0; // Local sum of lithiation contributions

// //     // Iterate over mesh elements efficiently
// //     for (int ei = 0; ei < nE; ++ei) {
// //         TmpF->GetNodalValues(ei, VtxVal); // Retrieve nodal values for the element
// //         double val = std::accumulate(VtxVal.begin(), VtxVal.end(), 0.0);
// //         EAvg(ei) = val / nC;
// //         lSum += EAvg(ei) * EVol(ei);
// //     }

// //     double gSum; 
// //     MPI_Allreduce(&lSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
// //     Xfr = gSum / gtps; // Calculate the degree of lithiation as the normalized sum
// // }

// // void Concentrations::CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value) {

// //     Rx2 = Rx1; // Copy the input reaction field to the output reaction field
// //     Rx2 *= value; // Scale the output reaction field by the specified factor

// // }

// // void Concentrations::CreateReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &Rx3, double value) {

// //     Rx3 = Rx1; // Copy the input reaction field to the output reaction field
// //     Rx3 += Rx2; // Add the second reaction field to the output reaction field
// //     Rx3 *= value; // Scale the output reaction field by the specified factor

// // }


// // void Concentrations::TotalReaction(mfem::ParGridFunction &Rx, double &xCrnt) {
    
// //     xCrnt = 0.0; // Initialize the local total reaction value to zero

// //     // calculate the west boundary size
// //     mfem::Vector Rmin, Rmax;
// //     gmesh->GetBoundingBox(Rmin, Rmax);

// //     int dim = gmesh->Dimension();
   
// //     if (dim == 2) L_w =  Rmax(1) - Rmin(1);
// //     else if (dim == 3) L_w = (Rmax(1) - Rmin(1))*(Rmax(2) - Rmin(2));

// //     for (int ei = 0; ei < nE; ei++) {
// //         Rx.GetNodalValues(ei, VtxVal); // Retrieve the nodal values of the reaction field for the current element
// //         double val = 0.0;
// //         for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}
// //         EAvg(ei) = val / nC; // Calculate the average reaction value for the element
// //         xCrnt += EAvg(ei) * EVol(ei); // Weight by the element volume and add to the total reaction
// //     }
    
// //     // Perform a global reduction to sum up contributions across all MPI processes
// //     MPI_Allreduce(&xCrnt, &geCrnt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
// //     // Calculate the reaction current density by normalizing with the domain characteristic length
// //     infx = geCrnt / (L_w);

// // }


// // void Concentrations::SaltConservation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx) {

// //     CeC = 0.0; // Initialize total salt concentration to zero
    
// //     // Temporary grid function to hold the product of concentration and potential fields
// //    static mfem::ParGridFunction CeT(fespace.get());

// //     // Calculate the product of concentration and potential for each element
// //     CeT = Cn;
// //     CeT *= psx;
    
// //     // Loop through all elements to calculate the average concentration for each element
// //     for (int ei = 0; ei < nE; ei++){
// //         CeT.GetNodalValues(ei,VtxVal); // Retrieve the nodal values for the current element 
// //         double val = 0.0;
// //         for (int vt = 0; vt < nC; vt++){val += VtxVal[vt];}        
// //         EAvg(ei) = val/nC; // Compute the average concentration for the element
// //         // Update the total salt concentration by adding the weighted average for the element
// //         CeC += EAvg(ei)*EVol(ei); 
// //     }
    
// //     // Perform a global reduction to sum up the salt concentration across all MPI processes
// //     MPI_Allreduce(&CeC, &gCeC, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);			
    
// //     // Compute the average concentration across the entire electrolyte
// //     CeAvg = gCeC/gtPse;	

// //     // Adjust the concentration field by normalizing the average concentration with the initial value
// //     Cn -= (CeAvg-Ce0);
// //     MPI_Barrier(MPI_COMM_WORLD);

// // }

