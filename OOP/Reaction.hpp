#ifndef REACTION_HPP
#define REACTION_HPP

#include "Concentrations_Base.hpp"
/**
 * @class Reaction
 * @brief Class for handling electrochemical reactions in battery simulations.
 *
 * This class provides methods to compute reaction rates, Butler-Volmer kinetics, 
 * exchange current density, and the total reaction current in a battery system.
 */
class Reaction {

public:

    /**
     * @brief Constructor for the Reaction class
     * @param pmesh Pointer to the parallel mesh
     * @param fespace Pointer to the finite element space
     * @param mh Reference to the mesh handler
     */
    Reaction(Initialize_Geometry &geo, Domain_Parameters &para);

    // virtual ~Reaction();

    Initialize_Geometry &geometry;
    Domain_Parameters &domain_parameters;

    /**
     * @brief Initializes the reaction field with a given initial value
     * @param Rx Reaction rate grid function
     * @param initial_value Initial value for the reaction rate
     */
    void Initialize(mfem::ParGridFunction &Rx, double initial_value);

    /**
     * @brief Computes the exchange current density and rate constants at the interface
     * @param Cn Concentration field at the interface
     */
    void ExchangeCurrentDensity(mfem::ParGridFunction &Cn);

    /**
     * @brief Applies the Butler-Volmer equation to calculate reaction rates
     * @param Rx Reaction rate grid function
     * @param Cn1 Concentration field on one side of the interface
     * @param Cn2 Concentration field on the other side of the interface
     * @param phx1 Potential field on one side of the interface
     * @param phx2 Potential field on the other side of the interface
     */
    void ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2);

    /**
     * @brief Computes the total reaction current over the domain
     * @param Rx Reaction rate grid function
     * @param global_current Output: Global reaction current across the domain
     */
    void TotalReactionCurrent(mfem::ParGridFunction &Rx, double &global_current);

    double global_current;     ///< Global reaction current



private:

    // MeshHandler &mesh_handler; ///< Reference to the mesh handler for geometry data


    /**
     * @brief Sets an initial reaction rate for the reaction field
     * @param Cn Concentration field
     * @param initial_value Initial reaction value
     */
    void SetInitialReaction(mfem::ParGridFunction &Cn, double initial_value);



    mfem::ParMesh *pmesh; ///< Pointer to the parallel mesh
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Pointer to the finite element space

    mfem::ParGridFunction AvP; ///< Grid function for active particle surface area
    mfem::ParGridFunction AvB; ///< Grid function for active boundary area

    int nE;                                         ///< Number of elements in the mesh
    int nC;                                         ///< Number of corners per element
    int nV;                                         ///< Number of vertices in the mesh

    double local_current; ///< Local reaction current for each MPI process

    // mfem::ParGridFunction *i0C; ///< Exchange current density field
    // mfem::ParGridFunction *OCV; ///< Open circuit voltage field
    // mfem::ParGridFunction *Kfw; ///< Open circuit voltage field
    // mfem::ParGridFunction *Kbw; ///< Backward reaction rate constant field
    // mfem::ParGridFunction *dPHE; ///< Voltage drop field

    std::unique_ptr<mfem::ParGridFunction> i0C;
    std::unique_ptr<mfem::ParGridFunction> OCV;
    std::unique_ptr<mfem::ParGridFunction> Kfw;
    std::unique_ptr<mfem::ParGridFunction> Kbw;
    std::unique_ptr<mfem::ParGridFunction> dPHE;

    const mfem::Vector& EVol; ///< Element volumes from the mesh handler


};

#endif // REACTION_HPP