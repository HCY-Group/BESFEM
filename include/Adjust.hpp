#ifndef ADJUST_HPP
#define ADJUST_HPP

#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "PotA.hpp"
#include "PotC.hpp"
#include "PotE.hpp"

/**
 * @class Adjust
 * @brief Handles voltage and current control for BESFEM battery simulations.
 *
 * The Adjust class provides routines for controlling electrochemical
 * boundary conditions during time-dependent simulations. Specifically, it
 * updates the applied electrode potentials to match the target global current
 * or to follow a specified voltage scan rate
 *
 */
class Adjust {

public:
    /**
     * @brief Construct an Adjust object.
     *
     * @param geo  Reference to the initialized geometry handler (meshes, FE space).
     * @param para Reference to the simulation domain parameter structure.
     */
    Adjust(Initialize_Geometry &geo, Domain_Parameters &para);

    /// Reference to geometry information (parallel mesh, FE space, ψ-fields).
    Initialize_Geometry &geometry;

    /// Reference to global domain parameters (material constants, scan rates, etc.).
    Domain_Parameters &domain_parameters;


    /**
     * @brief Adjusts to satisfy target global current for constant current simulations.
     *
     * After computing the updated surface potentials, this function updates
     * the associated grid functions (`phA_gf`, `phC_gf`) used in the FEM solve.
     *
     * @param current_A  The computed anode current contribution (A).
     * @param current_C  The computed cathode current contribution (A).
     * @param anode_potential   Reference to the anode potential object (PotA).
     * @param cathode_potential Reference to the cathode potential object (PotC).
     * @param phA_gf The ParGridFunction storing the anodic potential field.
     * @param phC_gf The ParGridFunction storing the cathodic potential field.
     * @param VCell  Output: computed or updated cell voltage (V).
     */
    void AdjustConstantCurrent(double current_A, double current_C, PotA &anode_potential, PotC &cathode_potential,
                              mfem::ParGridFunction &phA_gf, mfem::ParGridFunction &phC_gf, double &VCell);

private:

    /// Pointer to the parallel mesh used for the potential fields.
    mfem::ParMesh *pmesh;

    /// Shared pointer to the parallel finite element space for potentials.
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;
};

#endif // ADJUST_HPP
