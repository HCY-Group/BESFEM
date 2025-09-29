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
     * @brief Computes the exchange current density and rate constants at the interface using interpolation tables
     * @param Cn Concentration field at the interface
     */
    void TableExchangeCurrentDensity(mfem::ParGridFunction &Cn);

    /**
     * @brief Computes the exchange current density and rate constants at the interface
     * @param Cn Concentration field at the interface
     */
    void ExchangeCurrentDensity(mfem::ParGridFunction &Cn);


    /**
     * @brief Computes the exchange current density and rate constants at the interface using two concentration fields
     * @param Cn1 Concentration field on one side of the interface
     * @param Cn2 Concentration field on the other side of the interface
     */
    void ExchangeCurrentDensity(mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2);



    /**
     * @brief Computes the Butler-Volmer reaction rates based on concentration and potential fields
     * @param Rx Reaction rate grid function
     * @param Cn1 Concentration field on one side of the interface
     * @param Cn2 Concentration field on the other side of the interface
     * @param phx1 Potential field on one side of the interface
     * @param phx2 Potential field on the other side of the interface
     */
    void ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2);
    
    /**
     * @brief Computes the Butler-Volmer reaction rates based on concentration and potential fields for both anode and cathode
     * @param Rx Reaction rate grid function
     * @param Rx1 Reaction rate grid function for cathode
     * @param Rx2 Reaction rate grid function for anode
     * @param Cn1 Concentration field on cathode side of the interface
     * @param Cn2 Concentration field on anode side of the interface
     * @param Cn3 Concentration field in the electrolyte
     * @param phx1 Potential field on cathode side of the interface
     * @param phx2 Potential field on anode side of the interface
     * @param phx3 Potential field in the electrolyte
     */
    void ButlerVolmer(mfem::ParGridFunction &Rx, mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &Cn1, mfem::ParGridFunction &Cn2, mfem::ParGridFunction &Cn3, mfem::ParGridFunction &phx1, mfem::ParGridFunction &phx2, mfem::ParGridFunction &phx3);
    
    
    /**
     * @brief Computes the total reaction current across the domain
     * @param Rx Reaction rate grid function
     * @param global_current Output: Global reaction current across the domain
     */
    void TotalReactionCurrent(mfem::ParGridFunction &Rx, double &global_current);

    /**
     * @brief Applies the Butler-Volmer equation to calculate reaction rates
     * @param Rx Reaction rate grid function
     * @param Cn1 Concentration field on one side of the interface
     * @param Cn2 Concentration field on the other side of the interface
     * @param phx1 Potential field on one side of the interface
     * @param phx2 Potential field on the other side of the interface
     */

    std::unique_ptr<mfem::ParGridFunction> Kfw; ///< Open circuit voltage field
    std::unique_ptr<mfem::ParGridFunction> Kbw; ///< Backward reaction rate constant field
    std::unique_ptr<mfem::ParGridFunction> KfA; ///< Open circuit voltage field (anode)
    std::unique_ptr<mfem::ParGridFunction> KbA; ///< Backward reaction rate constant field (anode)
    std::unique_ptr<mfem::ParGridFunction> KfC; ///< Open circuit voltage field (cathode)
    std::unique_ptr<mfem::ParGridFunction> KbC; ///< Backward reaction rate constant field (cathode)
    std::unique_ptr<mfem::ParGridFunction> i0C; ///< Exchange current density field
    std::unique_ptr<mfem::ParGridFunction> OCV; ///< Open circuit voltage field

    // non-owning
    const mfem::ParGridFunction* AvP = nullptr;
    const mfem::ParGridFunction* AvB = nullptr;
    const mfem::ParGridFunction* AvA = nullptr;
    const mfem::ParGridFunction* AvC = nullptr;


private:

    /**
     * @brief Sets an initial reaction rate for the reaction field
     * @param Cn Concentration field
     * @param initial_value Initial reaction value
     */
    void SetInitialReaction(mfem::ParGridFunction &Cn, double initial_value);

    mfem::ParMesh *pmesh; ///< Pointer to the parallel mesh
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Pointer to the finite element space


    int nE;                                         ///< Number of elements in the mesh
    int nC;                                         ///< Number of corners per element
    int nV;                                         ///< Number of vertices in the mesh

    double local_current; ///< Local reaction current for each MPI process
    std::unique_ptr<mfem::ParGridFunction> dPHE; ///< Voltage drop field
    std::unique_ptr<mfem::ParGridFunction> dPHA; ///< Voltage drop field
    std::unique_ptr<mfem::ParGridFunction> dPHC; ///< Voltage

    const mfem::Vector& EVol; ///< Element volumes from the mesh handler

    // Interpolation tables
    mfem::Vector Ticks = mfem::Vector(101); ///< Ticks for interpolation
    mfem::Vector chmPot = mfem::Vector(101); ///< Charge transfer potential values
    mfem::Vector Mobility = mfem::Vector(101); ///< Mobility values for the Butler-Volmer equation
    mfem::Vector OCV_file = mfem::Vector(101); ///< Open circuit voltage values from file
    mfem::Vector i0_file = mfem::Vector(101); ///< Exchange current density values from file

    /**
     * @brief Gets the value from a table based on the concentration and ticks.
     * @param cn Concentration value to look up.
     * @param ticks Vector of tick values for interpolation.
     * @param data Vector of data values corresponding to the ticks.
     * @return Interpolated value from the table.
     */
    double GetTableValues(double cn, const mfem::Vector &ticks, const mfem::Vector &data);



};

#endif // REACTION_HPP