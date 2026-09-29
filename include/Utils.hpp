/**
 * @file Utils.hpp
 * @brief Defines utility routines for BESFEM field initialization, diagnostics, and output.
 */

#pragma once

#include <string>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <memory>
#include <vector>

#include "mfem.hpp"
#include "Initialize_Geometry.hpp"
#include "Domain_Parameters.hpp"

struct SimulationState;

/**
 * @class Utils
 * @brief Helper class for common BESFEM operations.
 *
 * Utils provides shared routines for initializing fields, computing lithiation
 * and reaction currents, evaluating global errors, computing pairwise flux
 * terms, and saving simulation snapshots.
 */
class Utils
{
public:
    /**
     * @brief Construct a Utils helper object.
     *
     * @param geo Reference to the geometry handler.
     * @param para Reference to the domain-parameter object.
     * @param cfg Reference to the simulation configuration.
     */
    Utils(Initialize_Geometry &geo, Domain_Parameters &para, const SimulationConfig &cfg);

    /**
     * @brief Set every local field entry to a uniform initial value.
     *
     * @param[out] Cn Field to initialize (for example, concentration or potential).
     * @param initial_value Value assigned to the field.
     */
    void SetInitialValue(mfem::ParGridFunction &Cn, double initial_value);

    /**
     * @brief Initialize two reaction fields.
     *
     * @param[out] Rx1 First reaction field.
     * @param[out] Rx2 Second reaction field.
     * @param value Initial value.
     */
    void InitializeReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value);

    /**
     * @brief Initialize three reaction fields.
     *
     * @param[out] Rx1 First reaction field.
     * @param[out] Rx2 Second reaction field.
     * @param[out] Rx3 Third reaction field.
     * @param value Initial value.
     */
    void InitializeReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &Rx3, double value);

    /**
     * @brief Compute and store the phase-weighted global concentration fraction.
     *
     * Integrates Cn * psx using element nodal averages and element volumes,
     * sums across MPI ranks, and divides by gtps. Retrieve with GetLithiation().
     * @pre gtps is nonzero; all ranks participate in the MPI reduction.
     *
     * @param Cn Concentration field.
     * @param psx Phase-field mask.
     * @param gtps Global integral of the phase-field mask.
     */
    void CalculateLithiation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, double gtps);

    /**
     * @brief Integrate the reaction field and normalize by transverse mesh extent.
     *
     * Sums element-average reaction values times element volumes across MPI
     * ranks, then divides by the bounding-box y extent (y times z in 3D).
     * @pre The transverse extent is nonzero; all ranks participate.
     *
     * @param Rx Reaction field.
     * @param[out] xCrnt Normalized global reaction integral.
     */
    void CalculateReactionInfx(mfem::ParGridFunction &Rx, double &xCrnt);

    /**
     * @brief Compute pairwise particle flux from chemical-potential differences.
     *
     * At each vertex, assigns weight * grad_psi * rho * Constants::Perm /
     * Constants::RT * (mu_nbr - mu_self).
     *
     * @param[out] sum_part Pair flux field, overwritten rather than accumulated.
     * @param weight Pairwise coupling weight.
     * @param grad_psi Pairwise interface/gradient field.
     * @param mu_self Chemical potential of the current particle.
     * @param mu_nbr Chemical potential of the neighboring particle.
     * @param rho Site density.
     */
    void ComputePairFlux(mfem::ParGridFunction &sum_part, mfem::ParGridFunction &weight, mfem::ParGridFunction &grad_psi,
                         mfem::ParGridFunction &mu_self, mfem::ParGridFunction &mu_nbr, double rho);

    /**
     * @brief Compute the normalized global sum of phase-weighted squared changes.
     *
     * Sums (px0 - potential)^2 * psx entries across MPI ranks and divides by
     * gtPsx. No square root or element-volume weighting is applied.
     * @pre gtPsx is nonzero; all ranks participate in the MPI reduction.
     *
     * @param px0 Previous potential field.
     * @param potential Updated potential field.
     * @param psx Phase-field mask.
     * @param[out] globalerror Normalized squared-change diagnostic.
     * @param gtPsx Global integral of the phase-field mask.
     */
    void CalculateGlobalError(mfem::ParGridFunction &px0, mfem::ParGridFunction &potential,
                              mfem::ParGridFunction &psx, double &globalerror, double gtPsx);

    /**
     * @brief Return the most recently computed lithiation fraction.
     *
     * @return Lithiation or concentration fraction.
     */
    double GetLithiation() const { return Xfr_; }

    /**
     * @brief Build a timestamped output directory path using local time.
     *
     * Returns ../outputs/Results/YYYYMMDD_HHMMSS without creating directories.
     *
     * @return Output directory path.
     */
    static inline std::string BuildRunOutdir()
    {

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

        std::ostringstream od;
        od << "../outputs/Results/" << ts.str();

        return od.str();
    }

    /**
     * @brief Save half-cell fields on timesteps divisible by save_interval.
     *
     * Writes concentration, potential, and combined reaction fields. At t = 0,
     * also writes the mesh and phase masks. Uses parallel SaveAsOne operations.
     * @param t Current timestep index.
     * @param outdir Existing output directory for snapshot files.
     * @param geometry Mesh and finite-element space.
     * @param domain_parameters Half-cell phase masks and particle-group masks.
     * @param state Initialized half-cell fields and electrode state.
     * @param electrode Active electrode (ANODE or CATHODE).
     * @param save_interval Number of timesteps between snapshots.
     * @pre save_interval is positive; all participating ranks call this method.
     */
    static void SaveHalfCellSnapshot(int t, const std::string& outdir, Initialize_Geometry& geometry,
        Domain_Parameters& domain_parameters, SimulationState& state, sim::Electrode electrode, int save_interval);

    /**
     * @brief Save full-cell fields on timesteps divisible by save_interval.
     *
     * Writes combined particle concentrations, electrolyte concentration,
     * potentials, and reaction fields. At t = 0, also writes the mesh and masks.
     * Uses parallel SaveAsOne operations.
     * @param t Current timestep index.
     * @param outdir Existing output directory for snapshot files.
     * @param geometry Mesh and finite-element space.
     * @param domain_parameters Anode, cathode, electrolyte, and particle masks.
     * @param state Initialized full-cell fields and electrode states.
     * @param save_interval Number of timesteps between snapshots.
     * @pre save_interval is positive; all participating ranks call this method.
     */
    static void SaveFullCellSnapshot(int t, const std::string& outdir, Initialize_Geometry& geometry,
        Domain_Parameters& domain_parameters, SimulationState& state, int save_interval);

    /**
     * @brief Print selected simulation settings and the output path on rank zero.
     * @param cfg Configuration supplying dt, dh, gc, Cr, and Vsr0.
     * @param outdir Output directory path to display.
     */
    static void PrintSimulationParameters(const SimulationConfig &cfg, const std::string &outdir);

    /**
     * @brief Print half-cell currents, field diagnostics, and lithiation on rank zero.
     * @param t Current timestep index.
     * @param VCell Cell voltage to report.
     * @param total_current Total current for the active electrode.
     * @param total_target Target current for the active electrode.
     * @param particle_currents Currents in active-electrode particle-group order.
     * @param state Initialized half-cell state supplying fields and lithiation.
     * @param para Domain data supplying group targets and phase totals.
     * @param electrode Active electrode (ANODE or CATHODE).
     * @pre particle_currents contains one entry per active particle group.
     */
    static void PrintHalfCellStatus(int t, double VCell, double total_current, double total_target,
        const std::vector<double> &particle_currents, const SimulationState &state, const Domain_Parameters &para, sim::Electrode electrode);

    /**
     * @brief Print full-cell currents, voltages, and lithiation on rank zero.
     *
     * Electrode-average lithiation is weighted by each group's global phase total.
     * @param t Current timestep index.
     * @param VCell Cell voltage to report.
     * @param anode_current Total anode current.
     * @param cathode_current Total cathode current.
     * @param state Initialized full-cell state supplying voltages and lithiation.
     * @param para Domain data supplying phase totals and target current.
     */
    static void PrintFullCellStatus(int t, double VCell, double anode_current, double cathode_current,
        const SimulationState &state, const Domain_Parameters &para);

    /**
     * @brief Print elapsed program time in whole seconds on rank zero.
     * @param start Timestamp recorded at program start.
     * @param end Timestamp recorded at program completion.
     */
    static void PrintProgramTime(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end);

    /**
     * @brief Check the configured timestep or voltage stopping condition.
     *
     * STEPS mode stops at t >= num_timesteps. VOLTAGE mode stops at
     * VCell <= VCut for positive Cr or VCell >= VCut for negative Cr.
     * A zero Cr does not trigger a voltage stop.
     * @param cfg Configuration specifying stop mode, step limit, Cr, and VCut.
     * @param t Current timestep index.
     * @param VCell Current cell voltage.
     * @return True when the selected stopping condition is met; false otherwise.
     */
    static bool ShouldStopSimulation(const SimulationConfig& cfg, int t, double VCell);

private:
    Initialize_Geometry &geometry_; ///< Geometry handler.
    Domain_Parameters   &domain_;   ///< Domain-parameter object.

    const SimulationConfig &cfg; ///< Simulation configuration.

    mfem::ParMesh *pmesh_ = nullptr; ///< Parallel mesh pointer.
    std::shared_ptr<mfem::ParFiniteElementSpace> fes_; ///< Parallel finite element space.

    int nE_ = 0; ///< Number of elements.
    int nC_ = 0; ///< Number of nodes per element.
    int nV_ = 0; ///< Number of vertices.

    mfem::Vector EVol_; ///< Element volumes.
    mfem::Vector EAvg_; ///< Per-element average workspace.
    mfem::Array<double> VtxVal_; ///< Vertex-value workspace.

    mfem::ParGridFunction TmpF_; ///< Temporary grid-function workspace.

    double Xfr_ = 0.0;   ///< Most recently computed lithiation fraction.
    double geCrnt_ = 0.0; ///< Global reaction current.
    double infx_ = 0.0;   ///< Integrated flux/current diagnostic.
};
