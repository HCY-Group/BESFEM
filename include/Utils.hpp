#pragma once
#include <string>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <ctime>

#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"

/**
 * @file Utils.hpp
 * @brief Utility helper class for BESFEM simulations.
 *
 * Provides common operations used across concentration and potential solvers:
 * initialization helpers, reduction-based metrics, lithiation tracking,
 * global error evaluation, and utilities for saving simulation snapshots.
 *
 * Also includes static helpers for constructing output directories based on
 * timestamps and mesh names.
 */

/**
 * @class Utils
 * @brief Convenience utilities for reaction, concentration, and potential updates.
 *
 * Responsibilities include:
 * - Initializing fields (concentration, reaction, potential)
 * - Computing lithiation (Xfr) and reaction currents
 * - Global MPI reductions and error metrics
 * - Saving simulation states to disk
 * - Building timestamped output directories
 *
 * Instances of Utils hold references to geometry and domain parameters needed
 * for consistent integration and reduction operations.
 */
class Utils
{
public:

    /**
     * @brief Construct a Utils helper instance.
     *
     * Stores references to geometry and domain parameters and initializes
     * local mesh metadata used by various utility routines.
     *
     * @param geo  Geometry handler (mesh, FE spaces).
     * @param para Domain parameter container.
     */
    Utils(Initialize_Geometry &geo, Domain_Parameters &para);

    // ----------------------------------------------------------------------
    // High-level field initialization & reduction utilities
    // ----------------------------------------------------------------------

    /**
     * @brief Assign a uniform initial value to a concentration field.
     *
     * @param Cn            Concentration grid function (in/out).
     * @param initial_value Value to assign everywhere in the domain.
     */
    void SetInitialValue(mfem::ParGridFunction &Cn, double initial_value);

    /**
     * @brief Initialize two reaction fields to a uniform value (half cell).
     *
     * @param Rx1 First reaction field.
     * @param Rx2 Second reaction field.
     * @param value Initial scalar value.
     */
    void InitializeReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value);

    /**
     * @brief Initialize three reaction fields to a uniform value (full cell).
     *
     * Useful in full-cell simulations where anode, cathode, and electrolyte
     * each carry their own reaction contributions.
     *
     * @param Rx1 Reaction field 1.
     * @param Rx2 Reaction field 2.
     * @param Rx3 Reaction field 3.
     * @param value Initial scalar value.
     */
    void InitializeReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &Rx3, double value);

    /**
     * @brief Compute the lithiation fraction.
     *
     * @param Cn   Concentration field.
     * @param psx  Phase mask ψ.
     * @param gtps Global integral of ψ (denominator).
     */
    void CalculateLithiation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, double gtps);

    /**
     * @brief Compute the total reaction current across the domain.
     *
     *
     * @param Rx     Reaction field.
     * @param xCrnt  Output: global reaction current.
     */
    void CalculateReactionInfx(mfem::ParGridFunction &Rx, double &xCrnt);

    /**
     * @brief Compute the global L2/RMS error between two potential fields.
     *
     * The error is masked by ψ and normalized by its global integral.
     *
     * @param px0        Previous potential field.
     * @param potential  Updated potential field.
     * @param psx        Phase mask ψ.
     * @param globalerror Output: global error (MPI reduced).
     * @param gtPsx       Global integral of ψ.
     */
    void CalculateGlobalError(mfem::ParGridFunction &px0, mfem::ParGridFunction &potential,
                              mfem::ParGridFunction &psx, double &globalerror, double gtPsx);

    /**
     * @brief Return the most recently computed lithiation fraction.
     */
    double GetLithiation() const { return Xfr_; }

    // ----------------------------------------------------------------------
    // Output directory builder
    // ----------------------------------------------------------------------

    /**
     * @brief Construct a timestamped output directory name.
     *
     * Format:
     * `../outputs/Results/YYYYMMDD_HHMMSS__nsteps=N__mesh=MeshName`
     *
     * @param mesh_file Path to the mesh file.
     * @param num_steps Number of time steps.
     * @return A formatted directory string.
     */
    static inline std::string BuildRunOutdir(const char* mesh_file, int num_steps)
    {
        namespace fs = std::filesystem;

        auto now = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);
        std::tm tm{};

    #if defined(_WIN32)
        localtime_s(&tm, &now_c);
    #else
        localtime_r(&now_c, &tm);
    #endif

        std::ostringstream ts;
        ts << std::put_time(&tm, "%Y%m%d_%H%M%S");
        std::string mesh_name = fs::path(mesh_file).stem().string();

        std::ostringstream od;
        od << "../outputs/Results/" << ts.str()
           << "__nsteps=" << num_steps
           << "__mesh=" << mesh_name;

        return od.str();
    }

    // ----------------------------------------------------------------------
    // Snapshot utilities (field dumps)
    // ----------------------------------------------------------------------

    /**
     * @brief Save a full simulation state snapshot (full-cell version).
     *
     * Dumps concentrations, potentials, ψ-weighted concentrations, and
     * particle-phase concentration CnP at a given timestep.
     *
     * @param t     Timestep index.
     * @param outdir Output directory.
     * @param geometry Initialized geometry container.
     * @param domain_parameters Domain parameters.
     * @param phA   Anode potential.
     * @param phC   Cathode potential.
     * @param phE   Electrolyte potential.
     * @param CnA   Anode concentration.
     * @param CnC   Cathode concentration.
     * @param CnE   Electrolyte concentration.
     * @param CnApsi ψ-weighted anode concentration.
     * @param CnCpsi ψ-weighted cathode concentration.
     * @param CnEpsi ψ-weighted electrolyte concentration.
     * @param CnP   Particle concentration (if present).
     * @param save_interval Frequency of saving.
     */
    static void SaveSimulationSnapshot(int t, const std::string &outdir,
        Initialize_Geometry &geometry, Domain_Parameters &domain_parameters,
        mfem::ParGridFunction &phA, mfem::ParGridFunction &phC, mfem::ParGridFunction &phE,
        mfem::ParGridFunction &CnA, mfem::ParGridFunction &CnC, mfem::ParGridFunction &CnE,
        mfem::ParGridFunction &CnApsi, mfem::ParGridFunction &CnCpsi, mfem::ParGridFunction &CnEpsi,
        mfem::ParGridFunction &CnP, int save_interval = 500);

    /**
     * @brief Save a reduced simulation snapshot (half-cell version).
     *
     * Only cathode + electrolyte fields are stored.
     *
     * @param t     Timestep index.
     * @param outdir Output directory.
     * @param geometry Initialized geometry container.
     * @param domain_parameters Domain parameters.
     * @param phC   Cathode potential.
     * @param phE   Electrolyte potential.
     * @param CnC   Cathode concentration.
     * @param CnE   Electrolyte concentration.
     * @param CnCpsi ψ-weighted cathode concentration.
     * @param CnEpsi ψ-weighted electrolyte concentration.
     * @param save_interval Frequency of saving.
     */
    static void SaveSimulationSnapshot(int t, const std::string &outdir,
        Initialize_Geometry &geometry, Domain_Parameters &domain_parameters,
        mfem::ParGridFunction &phC, mfem::ParGridFunction &phE,
        mfem::ParGridFunction &CnC, mfem::ParGridFunction &CnE,
        mfem::ParGridFunction &CnCpsi, mfem::ParGridFunction &CnEpsi,
        int save_interval = 500);

    static void SaveSimulationSnapshotMulti(int t, const std::string &outdir,
    Initialize_Geometry &geometry, Domain_Parameters &domain_parameters,
    mfem::ParGridFunction &CnC_1, mfem::ParGridFunction &CnC_2, mfem::ParGridFunction &CnC_3,
    mfem::ParGridFunction &C1_out, mfem::ParGridFunction &C2_out, mfem::ParGridFunction &C3_out,
    int save_interval = 500);

private:

    // Stored references
    Initialize_Geometry &geometry_;
    Domain_Parameters   &domain_;

    // Cached mesh data
    mfem::ParMesh *pmesh_;
    std::shared_ptr<mfem::ParFiniteElementSpace> fes_;

    int nE_, nC_, nV_;
    mfem::Vector EVol_;
    mfem::Vector EAvg_;
    mfem::Array<double> VtxVal_;

    mfem::ParGridFunction TmpF_;

    // Cached global quantities
    double Xfr_   = 0.0; ///< Lithiation fraction.
    double geCrnt_ = 0.0; ///< Global reaction current.
    double infx_   = 0.0; ///< Boundary flux.
};
