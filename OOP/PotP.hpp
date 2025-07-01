#ifndef POTP_HPP
#define POTP_HPP

#include "Potentials_Base.hpp"

/**
 * @brief Boundary value for the particle potential
 */
extern double BvP;

/**
 * @class PotP
 * @brief Derived class implementing the particle potential model for battery simulations
 *
 * This class extends the `Potentials` base class to handle computations
 * related to particle potentials, including initialization, time-stepping,
 * and global error calculations
 */
class PotP : public Potentials {

    // std::unique_ptr<mfem::ParBilinearForm> K; ///< Unique pointer for the bilinear form used in matrix assembly



public:

    // /**
    //  * @brief Constructor for the PotP class
    //  * @param pm Pointer to the parallel mesh
    //  * @param fe Pointer to the finite element space
    //  * @param mh Reference to the mesh handler
    //  */
    PotP(Initialize_Geometry &geo, Domain_Parameters &para);

    // // virtual ~PotP();


    Initialize_Geometry &geometry;
    Domain_Parameters &domain_parameters;

    double BvP; ///< Boundary value for particle potential


    // /**
    //  * @brief Initializes the particle potential field and solver
    //  * @param ph Particle potential grid function
    //  * @param initial_value Initial potential value
    //  */
    void Initialize(mfem::ParGridFunction &Cn, double initial_value, mfem::ParGridFunction &psx);

    // /**
    //  * @brief Performs time-stepping for the particle potential
    //  * @param Cn Particle concentration field
    //  * @param psx Psi potential field
    //  * @param phx Particle potential field
    //  */
    void TimeStep(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, mfem::ParGridFunction &potential);
    void Advance(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);

    // /**
    //  * @brief Calculates the global error in the particle potential solution
    //  * @param Rx Reaction field
    //  * @param phx Particle potential field
    //  * @param psx Psi potential field
    //  * @param gerror Global error (output)
    //  */
    // void CalculateGlobalError(mfem::ParGridFunction &Rx, mfem::ParGridFunction &phx, mfem::ParGridFunction &psx, double &gerror);


    // // static mfem::CGSolver *cgPP_solver; ///< Static variable for the conjugate gradient solver
    // std::unique_ptr<mfem::CGSolver> cgPP_solver;

    // // /**
    // //  * @brief Getter for the static conjugate gradient solver instance
    // //  * @return Pointer to the conjugate gradient solver
    // //  */
    // // static mfem::CGSolver *GetcgPPsolver() { return cgPP_solver; }

    // double error_P = 1.0; ///< Local error in the particle potential solution



private:

    void ParticleConductivity(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx); ///< Computes particle conductivity

    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Pointer to the finite element space

    std::shared_ptr<mfem::CGSolver> cgPP_solver;

    mfem::ParLinearForm B1t; ///< Linear form for the right-hand side
    mfem::HypreParVector X1v; ///< Solution vector
    mfem::HypreParVector B1v; ///< Right-hand-side vector
    mfem::HypreParVector Fpb; ///< Right-hand-side vector
    mfem::HypreParVector Xs0; ///< Solution vector for particle potential
    mfem::HypreSmoother Mpp; ///< Preconditioner for the solver

    double gtPsi; ///< Total Psi from MeshHandler

    mfem::Array<int> ess_tdof_list_e; ///< List of essential true degrees of freedom for Dirichlet boundary conditions
    mfem::Array<int> dbc_e_bdr; ///< Array marking Dirichlet boundary attributes

    mfem::ParGridFunction kap; ///< Conductivity field for particle potential
    mfem::GridFunctionCoefficient cKp; ///< Coefficient for the conductivity field
    mfem::GridFunctionCoefficient cRp; ///< Coefficient for the reaction field
    mfem::ParGridFunction RpP; ///< Reaction field for the particle potential
    mfem::ParGridFunction pP0; ///< Particle potential grid function

    std::unique_ptr<mfem::ParBilinearForm> Kp2; ///< Bilinear form for particle potential conductivity
    std::shared_ptr<mfem::HypreParMatrix> KmP; ///< Stiffness matrix for particle potential conductivity
    std::unique_ptr<mfem::ParLinearForm> Bp2; ///< Linear form for the reaction term

    mfem::ParLinearForm Fpt; ///< Linear form for the force term in particle potential calculations

};

#endif // POTP_HPP