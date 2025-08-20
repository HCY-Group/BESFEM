#ifndef CURRENT_HPP
#define CURRENT_HPP


#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"

/**
 * @class Current
 * @brief Class for handling current control in battery simulations
 *
 * This class provides methods to adjust the particle potential based on the target
 * global current and voltage scan rate during simulations
 */
class Current {

public:
    /**
     * @brief Constructor for the Current class
     *
     * Initializes the Current object and associates it with the parallel mesh,
     * finite element space, and mesh handler
     *
     * @param pmesh Pointer to the parallel mesh
     * @param fespace Pointer to the finite element space
     * @param mh Reference to the mesh handler
     */
    Current(Initialize_Geometry &geo, Domain_Parameters &para);

    Initialize_Geometry &geometry;
    Domain_Parameters &domain_parameters;
    
    /**
     * @brief Adjusts the particle potential to maintain the constant target global current
     *
     * Updates the particle potential and adjusts the boundary value of the particle
     * potential (`BvP`) based on the voltage scan rate (`Vsr`) and the time step (`dt`)
     *
     * @param phx Particle potential grid function
     * @param global_current Current value computed across the entire domain
     */
    void Constant(mfem::ParGridFunction &phx, double &global_current);


private:

    mfem::ParMesh *pmesh; ///< Pointer to the parallel mesh
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace; ///< Pointer to the finite element space

    double sgn; ///< Sign of the difference between the target and global current
    double dV; ///< Voltage adjustment step based on the target current
 
};


#endif // CURRENT_HPP