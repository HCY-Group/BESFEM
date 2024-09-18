// Reaction.hpp
#ifndef REACTION_HPP
#define REACTION_HPP

#include "mfem.hpp"
#include "Mesh_Handler.hpp"
#include "Constants.hpp"
#include "Concentrations.hpp"

#include <memory>

using namespace mfem;

class Concentrations; // Forward declaration

class Reaction {
public:
    Reaction(MeshHandler &mesh_handler, Concentrations &concentrations); // Constructor

    void Initialize(); // Method to initialize Rxn

    std::unique_ptr<mfem::ParGridFunction> Rxn;     // Rxn (ParGridFunction)

private:
    
    MeshHandler &mesh_handler;
    Concentrations &concentrations; // Reference to Concentrations
    // ParGridFunction Rxn; // ParGridFunction for Rxn

    void CreateRx(mfem::ParGridFunction &Rx, double initial_value);
    void SetAvP(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Av, double value);

    mfem::ParFiniteElementSpace* fespace;           // Finite element space

    // std::unique_ptr<mfem::ParGridFunction> Rxn;     // Rxn (ParGridFunction)



};

#endif // REACTION_HPP
