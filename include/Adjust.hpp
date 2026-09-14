#ifndef ADJUST_HPP
#define ADJUST_HPP

#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"
#include "ElectrodePotential.hpp"
#include "ElectrolytePotential.hpp"

class ElectrolytePotential;

/**
 * @class Adjust
 * @brief Controls applied voltage/current conditions during BESFEM simulations.
 *
 * The Adjust class updates electrode boundary potentials during time-dependent
 * battery simulations. It is mainly used for constant-current control, where
 * the applied anode/cathode potentials are adjusted so the simulated total
 * current approaches the target current.
 */
class Adjust {

public:
    /**
     * @brief Construct an Adjust object.
     *
     * Stores references to the geometry, domain parameters, and simulation
     * configuration needed for voltage/current control.
     *
     * @param geo Reference to the initialized geometry object.
     * @param para Reference to the domain parameter object.
     * @param cfg Reference to the simulation configuration.
     */
    Adjust(Initialize_Geometry &geo, Domain_Parameters &para, const SimulationConfig &cfg);

    /// Reference to geometry data, including the mesh and finite element space.
    Initialize_Geometry &geometry;

    /// Reference to global simulation fields and electrochemical parameters.
    Domain_Parameters &domain_parameters;

    /// Reference to user-defined simulation settings.
    const SimulationConfig &cfg;

    /**
     * @brief Adjust electrode potentials to satisfy a target global current.
     *
     * Compares the computed anode and cathode currents against the target
     * current stored in the domain parameters, then updates the applied
     * electrode potentials. The corresponding potential grid functions are
     * also updated so the next FEM solve uses the adjusted boundary values.
     *
     * @param current_A Computed anode current contribution.
     * @param current_C Computed cathode current contribution.
     * @param anode_potential Anode electrode-potential object.
     * @param cathode_potential Cathode electrode-potential object.
     * @param phA_gf Parallel grid function storing the anode solid potential.
     * @param phC_gf Parallel grid function storing the cathode solid potential.
     * @param VCell Output cell voltage after the adjustment.
     */
    void AdjustConstantCurrent(double current_A,
                               double current_C,
                               ElectrodePotential &anode_potential,
                               ElectrodePotential &cathode_potential,
                               mfem::ParGridFunction &phA_gf,
                               mfem::ParGridFunction &phC_gf,
                               double &VCell);


    void AdjustHalfCellCurrent(double total_current, double total_target, ElectrolytePotential& electrolyte_potential, mfem::ParGridFunction& phE);

private:
    /// Pointer to the parallel mesh used by the potential solve.
    mfem::ParMesh *pmesh;

    /// Parallel finite element space used by the potential fields.
    std::shared_ptr<mfem::ParFiniteElementSpace> fespace;
};

#endif // ADJUST_HPP