// /**
//  * @file Reaction.cpp
//  * @brief Implementation of the Reaction class for electrochemical reaction modeling in batteries.
//  */

// #include "Reaction.hpp"
// #include "../code/Constants.hpp"
// #include "mfem.hpp"
// #include "CnE.hpp"
// #include "PotE.hpp"
// #include "PotP.hpp"


// Reaction::Reaction(Initialize_Geometry &geo, Domain_Parameters &para)
//     : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para), AvP(*para.AvP), AvB(*para.AvB), EVol(para.EVol)

// {
    
//     nE = geometry.nE; 
//     nC = geometry.nC; 
//     nV = geometry.nV; 

//     i0C = std::make_unique<mfem::ParGridFunction>(fespace.get()); // exchange current density
//     OCV = std::make_unique<mfem::ParGridFunction>(fespace.get()); // open circuit voltage
//     Kfw = std::make_unique<mfem::ParGridFunction>(fespace.get()); // forward reaction constant
//     Kbw = std::make_unique<mfem::ParGridFunction>(fespace.get()); // backward rection constant
//     dPHE = std::make_unique<mfem::ParGridFunction>(fespace.get()); // voltage drop


// }

// double InterpolateFromTable(double cn, const mfem::Vector &X_101, const mfem::Vector &table_101)
// {
//     // Fortran-style clamping
//     if (cn < 0.0) cn = 1.0e-6;
//     if (cn > 1.0) cn = 1.0;

//     int idx = std::floor(cn / 0.01);
//     if (idx < 0) idx = 0;
//     if (idx > 99) idx = 99;

//     return table_101(idx) + (cn - X_101(idx)) / 0.01 * (table_101(idx + 1) - table_101(idx));
// }


// void Reaction::Initialize(mfem::ParGridFunction &Rx, double initial_value) {

//     SetInitialReaction(Rx, initial_value);
//     Rx = AvP; // Scale by active particle surface area
//     Rx *= 1.0e-8; // Apply a scaling factor

// }

// void Reaction::SetInitialReaction(mfem::ParGridFunction &Rx, double initial_value) {
//     for (int i = 0; i < Rx.Size(); ++i) {
//         Rx(i) = initial_value;
//     }
// }

// void Reaction::ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2){

//     for (int vi = 0; vi < nV; vi++){
//         if ( AvB(vi) * Constants::dh > 0.0 ){ // Check for interface presence

//             (*dPHE)(vi) = phx1(vi) - phx2(vi); // Voltage drop across the interface

//             Rx(vi) = AvP(vi) * ((*Kfw)(vi)*Cn2(vi)*exp(-Constants::alp*Constants::Cst1*(*dPHE)(vi)) - \
// 					                   (*Kbw)(vi)*Cn1(vi)*exp( Constants::alp*Constants::Cst1*(*dPHE)(vi)));


//         }
//     }
// }


// // rate constants and exchange current density at interface
// void Reaction::ExchangeCurrentDensity(mfem::ParGridFunction &Cn){
//     for (int vi = 0; vi < nV; vi++){
//         if(AvB(vi) * Constants::dh > 0.0){ 
//             double val = -0.2 * (Cn(vi) - 0.37) - 1.559 - 0.9376 * tanh(8.961 * Cn(vi) - 3.195);
//             (*i0C)(vi) = pow(10.0, val) * 1.0e-3; // Exchange current density
//             (*OCV)(vi) = 1.095 * Cn(vi) * Cn(vi) - 8.324e-7 * exp(14.31 * Cn(vi)) + 4.692 * exp(-0.5389 * Cn(vi)); // open circuit voltage
//             (*Kfw)(vi) = (*i0C)(vi) / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * (*OCV)(vi)); // forward reaction constant
//             (*Kbw)(vi) = (*i0C)(vi) / (Constants::Frd * Cn(vi)) * exp(-Constants::alp * Constants::Cst1 * (*OCV)(vi)); // backward rection constant
//         }
//     }
// }


// void Reaction::TotalReactionCurrent(mfem::ParGridFunction &Rx, double &global_current){

//     local_current = 0.0;
//     mfem::Array<double> VtxVal(nC);
//     mfem::Vector EAvg(nE);

//     // Loop over all elements to calculate local reaction current
//     for (int ei = 0; ei < nE; ei++){
//         Rx.GetNodalValues(ei,VtxVal) ;
//         // double val = 0.0;
//         // for (int vt = 0; vt < nC; vt++){
//         //     val += VtxVal[vt];
//         // }
//         double val = std::accumulate(VtxVal.begin(), VtxVal.end(), 0.0);
//         EAvg(ei) = val/nC;	// Average reaction rate for the element
//         local_current += EAvg(ei)*EVol(ei) ; // Accumulate local current		 			
//     }
    
//     // Perform global reduction to compute the total current
//     MPI_Allreduce(&local_current, &global_current, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

// }


// /**
//  * @file Reaction.cpp
//  * @brief Implementation of the Reaction class for electrochemical reaction modeling in batteries.
//  */

//  #include "Reaction.hpp"
//  #include "../code/Constants.hpp"
//  #include "mfem.hpp"
//  #include "CnE.hpp"
//  #include "PotE.hpp"
//  #include "PotP.hpp"
//  #include <fstream>
//  #include <cmath>
 
//  Reaction::Reaction(Initialize_Geometry &geo, Domain_Parameters &para)
//      : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo), domain_parameters(para), AvP(*para.AvP), AvB(*para.AvB), EVol(para.EVol)
//  {
//      nE = geometry.nE; 
//      nC = geometry.nC; 
//      nV = geometry.nV; 
 
//      i0C = std::make_unique<mfem::ParGridFunction>(fespace.get()); // exchange current density
//      OCV = std::make_unique<mfem::ParGridFunction>(fespace.get()); // open circuit voltage
//      Kfw = std::make_unique<mfem::ParGridFunction>(fespace.get()); // forward reaction constant
//      Kbw = std::make_unique<mfem::ParGridFunction>(fespace.get()); // backward rection constant
//      dPHE = std::make_unique<mfem::ParGridFunction>(fespace.get()); // voltage drop
 
//      // Load table data for interpolation
//      X_101.SetSize(101);
//      i0C_101.SetSize(101);
//      OCV_101.SetSize(101);
 
//      std::ifstream myXfile("../C_Li_X_101.txt");
//      std::ifstream myI0file("../C_Li_J2_101.txt");
//      std::ifstream myOCVfile("../C_Li_O3_101.txt");
 
//      for (int i = 0; i < 101; i++) myXfile >> X_101(i);
//      for (int i = 0; i < 101; i++) myI0file >> i0C_101(i);    
//      for (int i = 0; i < 101; i++) myOCVfile >> OCV_101(i);
//  }
 
//  double Reaction::InterpolateFromTable(double cn, const mfem::Vector &X_101, const mfem::Vector &table_101)
//  {
//      if (cn < 0.0) cn = 1.0e-6;
//      if (cn > 1.0) cn = 1.0;
 
//      int idx = std::floor(cn / 0.01);
//      if (idx < 0) idx = 0;
//      if (idx > 99) idx = 99;
 
//      return table_101(idx) + (cn - X_101(idx)) / 0.01 * (table_101(idx + 1) - table_101(idx));
//  }
 
//  void Reaction::Initialize(mfem::ParGridFunction &Rx, double initial_value) {
//      SetInitialReaction(Rx, initial_value);
//     //  Rx += AvP;
//     //  Rx *= 1.0e-8;
//  }
 
//  void Reaction::SetInitialReaction(mfem::ParGridFunction &Rx, double initial_value) {
//      for (int i = 0; i < Rx.Size(); ++i) {
//          Rx(i) = initial_value;
//      }
//  }
 
//  void Reaction::ButlerVolmer(mfem::ParGridFunction &Rx, const mfem::ParGridFunction &Cn1,
//                             const mfem::ParGridFunction &Cn2, const mfem::ParGridFunction &phx1,
//                             const mfem::ParGridFunction &phx2) {
//     for (int vi = 0; vi < nV; ++vi) {
//         if (AvP(vi) > 0.01) {
//             double eta = phx1(vi) - phx2(vi);
//             (*dPHE)(vi) = eta;

//             double kfw = (*Kfw)(vi);
//             double kbw = (*Kbw)(vi);
//             double cnp = std::max(Cn1(vi), 1.0e-6); // fix this
//             double cne = Cn2(vi);

//             double rxn = kfw * cne * std::exp(-Constants::alp * Constants::Cst1 * eta)
//             - kbw * cnp * std::exp(Constants::alp * Constants::Cst1 * eta);

//             Rx(vi) = rxn;
//             // Rx(vi) = AvP(vi) * rxn; // SBM: restrict reaction to interface

//         }
//     }
// }
    

 
//  void Reaction::ExchangeCurrentDensity(mfem::ParGridFunction &Cn){
//      for (int vi = 0; vi < nV; vi++){
//         //  if(AvB(vi) * Constants::dh > 0.0){
//              double cn_val = Cn(vi);
//              cn_val = std::max(cn_val, 1.0e-6);  // Avoid division by zero

//              double i0 = InterpolateFromTable(cn_val, X_101, i0C_101) * 1.0e-3; // Convert mA to A
//              double ocv = InterpolateFromTable(cn_val, X_101, OCV_101);
 
//              (*i0C)(vi) = i0;
//              (*OCV)(vi) = ocv;
//              (*Kfw)(vi) = i0 / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * ocv);
//              (*Kbw)(vi) = i0 / (Constants::Frd * cn_val) * exp(-Constants::alp * Constants::Cst1 * ocv);
//         //  }
//      }

//     //  for (int i = 0; i < Cn.Size(); i++) {
//     //     double cn_val = Cn(i);
//     //     (*i0C)(i) = InterpolateFromTable(cn_val, X_101, i0C_101) * 1.0e-3 ;
//     //     (*OCV)(i) = InterpolateFromTable(cn_val, X_101, OCV_101) ;
    
//     // }
    
//     // i0C->Save("Results2/ioc"); // save mobility values for debugging
//     // OCV->Save("Results2/ocv"); // save mobility values for debugging
//  }



 
//  void Reaction::TotalReactionCurrent(mfem::ParGridFunction &Rx, double &global_current){
//      local_current = 0.0;
//      mfem::Array<double> VtxVal(nC);
//      mfem::Vector EAvg(nE);
 
//      for (int ei = 0; ei < nE; ei++){
//          Rx.GetNodalValues(ei,VtxVal);
//          double val = std::accumulate(VtxVal.begin(), VtxVal.end(), 0.0);
//          EAvg(ei) = val/nC;
//          local_current += EAvg(ei)*EVol(ei);
//      }
 
//      MPI_Allreduce(&local_current, &global_current, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  }
 
// //  void Reaction::WriteKfwToFile(const std::string &filename) const {
// //     std::ofstream out(filename);
// //     if (!out.is_open()) {
// //         std::cerr << "Failed to open file " << filename << " for writing Kfw values.\n";
// //         return;
// //     }

// //     for (int i = 0; i < Kfw->Size(); ++i) {
// //         out << std::setprecision(12) << (*Kfw)(i) << "\n";
// //     }

// //     out.close();
// //     std::cout << "Kfw values written to " << filename << "\n";
// // }


/**
 * @file Reaction.cpp
 * @brief Implementation of the Reaction class for electrochemical reaction modeling in batteries.
 */

 #include "../include/Reaction.hpp"
 #include "../inputs/Constants.hpp"
 #include "mfem.hpp"
 #include "../include/CnE.hpp"
 #include "../include/PotE.hpp"
 #include "../include/PotP.hpp"
 #include <fstream>
 #include <cmath>
 
 Reaction::Reaction(Initialize_Geometry &geo, Domain_Parameters &para)
     : pmesh(geo.parallelMesh.get()), fespace(geo.parfespace), geometry(geo),
       domain_parameters(para), AvP(*para.AvP), AvB(*para.AvB), EVol(para.EVol)
 {
     nE = geometry.nE; 
     nC = geometry.nC; 
     nV = geometry.nV; 
 
     i0C = std::make_unique<mfem::ParGridFunction>(fespace.get()); // exchange current density
     OCV = std::make_unique<mfem::ParGridFunction>(fespace.get()); // open circuit voltage
     Kfw = std::make_unique<mfem::ParGridFunction>(fespace.get()); // forward reaction constant
     Kbw = std::make_unique<mfem::ParGridFunction>(fespace.get()); // backward reaction constant
     dPHE = std::make_unique<mfem::ParGridFunction>(fespace.get()); // voltage drop
 
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
    //  Rx = AvP; // Scale by active particle surface area
    //  Rx *= 8.0e-10 * 0.99 * 2.25; // Apply a scaling factor
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
        if(AvB(vi) * Constants::dh > 0.0){ 
            double val = -0.2 * (Cn(vi) - 0.37) - 1.559 - 0.9376 * tanh(8.961 * Cn(vi) - 3.195);
            (*i0C)(vi) = pow(10.0, val) * 1.0e-3; // Exchange current density
            (*OCV)(vi) = 1.095 * Cn(vi) * Cn(vi) - 8.324e-7 * exp(14.31 * Cn(vi)) + 4.692 * exp(-0.5389 * Cn(vi)); // open circuit voltage
            (*Kfw)(vi) = (*i0C)(vi) / (Constants::Frd * 0.001) * exp(Constants::alp * Constants::Cst1 * (*OCV)(vi)); // forward reaction constant
            (*Kbw)(vi) = (*i0C)(vi) / (Constants::Frd * Cn(vi)) * exp(-Constants::alp * Constants::Cst1 * (*OCV)(vi)); // backward rection constant
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
        if ( AvB(vi) * Constants::dh > 0.0 ){ // Check for interface presence
            (*dPHE)(vi) = phx1(vi) - phx2(vi); // Voltage drop across the interface
            Rx(vi) = AvP(vi) * ((*Kfw)(vi)*Cn2(vi)*exp(-Constants::alp*Constants::Cst1*(*dPHE)(vi)) - \
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
 