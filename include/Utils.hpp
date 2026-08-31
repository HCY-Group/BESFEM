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
     * @brief Set a concentration field to a uniform initial value.
     *
     * @param Cn Concentration field to initialize.
     * @param initial_value Value assigned to the field.
     */
    void SetInitialValue(mfem::ParGridFunction &Cn, double initial_value);

    /**
     * @brief Initialize two reaction fields.
     *
     * @param Rx1 First reaction field.
     * @param Rx2 Second reaction field.
     * @param value Initial value.
     */
    void InitializeReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value);

    /**
     * @brief Initialize three reaction fields.
     *
     * @param Rx1 First reaction field.
     * @param Rx2 Second reaction field.
     * @param Rx3 Third reaction field.
     * @param value Initial value.
     */
    void InitializeReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &Rx3, double value);

    /**
     * @brief Compute the global lithiation or concentration fraction.
     *
     * @param Cn Concentration field.
     * @param psx Phase-field mask.
     * @param gtps Global integral of the phase-field mask.
     */
    void CalculateLithiation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, double gtps);

    /**
     * @brief Compute the global reaction current.
     *
     * @param Rx Reaction field.
     * @param xCrnt Output global reaction current.
     */
    void CalculateReactionInfx(mfem::ParGridFunction &Rx, double &xCrnt);

    /**
     * @brief Compute pairwise particle flux from chemical-potential differences.
     *
     * @param sum_part Output accumulated pair flux field.
     * @param weight Pairwise coupling weight.
     * @param grad_psi Pairwise interface/gradient field.
     * @param mu_self Chemical potential of the current particle.
     * @param mu_nbr Chemical potential of the neighboring particle.
     * @param rho Site density.
     */
    void ComputePairFlux(mfem::ParGridFunction &sum_part, mfem::ParGridFunction &weight, mfem::ParGridFunction &grad_psi,
                         mfem::ParGridFunction &mu_self, mfem::ParGridFunction &mu_nbr, double rho);

    /**
     * @brief Compute the global error between two potential fields.
     *
     * @param px0 Previous potential field.
     * @param potential Updated potential field.
     * @param psx Phase-field mask.
     * @param globalerror Output global error.
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
     * @brief Build a timestamped output directory path.
     *
     * @return Output directory path.
     */
    static inline std::string BuildRunOutdir()
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

        std::ostringstream od;
        od << "../outputs/Results/" << ts.str();

        return od.str();
    }

    /**
     * @brief Save a full-cell simulation snapshot.
     *
     * @param t Timestep index.
     * @param outdir Output directory.
     * @param geometry Geometry handler.
     * @param domain_parameters Domain-parameter object.
     * @param phA Anode potential field.
     * @param phC Cathode potential field.
     * @param phE Electrolyte potential field.
     * @param CnA Anode concentration field.
     * @param CnC Cathode concentration field.
     * @param CnE Electrolyte concentration field.
     * @param CnApsi Masked anode concentration field.
     * @param CnCpsi Masked cathode concentration field.
     * @param CnEpsi Masked electrolyte concentration field.
     * @param CnP Combined particle concentration field.
     * @param save_interval Number of timesteps between saved snapshots.
     */
    static void SaveSimulationSnapshot(int t, const std::string &outdir,
        Initialize_Geometry &geometry, Domain_Parameters &domain_parameters,
        mfem::ParGridFunction &phA, mfem::ParGridFunction &phC, mfem::ParGridFunction &phE,
        mfem::ParGridFunction &CnA, mfem::ParGridFunction &CnC, mfem::ParGridFunction &CnE,
        mfem::ParGridFunction &CnApsi, mfem::ParGridFunction &CnCpsi, mfem::ParGridFunction &CnEpsi,
        mfem::ParGridFunction &CnP, int save_interval = 500);

    /**
     * @brief Save a half-cell simulation snapshot.
     *
     * @param t Timestep index.
     * @param outdir Output directory.
     * @param geometry Geometry handler.
     * @param domain_parameters Domain-parameter object.
     * @param phC Electrode potential field.
     * @param phE Electrolyte potential field.
     * @param CnC Electrode concentration field.
     * @param CnE Electrolyte concentration field.
     * @param CnCpsi Masked electrode concentration field.
     * @param CnEpsi Masked electrolyte concentration field.
     * @param save_interval Number of timesteps between saved snapshots.
     */
    static void SaveSimulationSnapshot(int t, const std::string &outdir,
        Initialize_Geometry &geometry, Domain_Parameters &domain_parameters,
        mfem::ParGridFunction &phC, mfem::ParGridFunction &phE,
        mfem::ParGridFunction &CnC, mfem::ParGridFunction &CnE,
        mfem::ParGridFunction &CnCpsi, mfem::ParGridFunction &CnEpsi,
        int save_interval = 500);


    /**
     * @brief Save multi-particle simulation snapshots.
     * 
     * @param t Timestep index.
     * @param outdir Output directory.
     * @param geometry Geometry handler.
     * @param domain_parameters Domain-parameter object.
     * @param particle_cn Particle concentration fields to save.
     * @param particle_ps Particle phase fields.
     * @param electrode_psi Electrode phase field.
     * @param particle_out Output workspaces for masked/saved fields.
     * @param electrode_name Name of the electrode.
     * @param save_interval Number of timesteps between saved snapshots.
     */
    static void SaveSimulationSnapshotMulti(int t, const std::string& outdir,
        Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
        const std::vector<mfem::ParGridFunction*>& particle_cn,
        const std::vector<std::unique_ptr<mfem::ParGridFunction>>& particle_ps,
        mfem::ParGridFunction& electrode_psi,
        std::vector<std::unique_ptr<mfem::ParGridFunction>>& particle_out,
        const std::string& electrode_name, int save_interval);

    static void SaveCombinedSnapshot(int t, const std::string& outdir,
        Initialize_Geometry& geometry, Domain_Parameters& domain_parameters,
        const std::vector<mfem::ParGridFunction*>& anode_fields,
        const std::vector<mfem::ParGridFunction*>& anode_ps,
        const std::vector<mfem::ParGridFunction*>& cathode_fields,
        const std::vector<mfem::ParGridFunction*>& cathode_ps,
        const std::vector<mfem::ParGridFunction*>& electrolyte_fields,
        const std::vector<mfem::ParGridFunction*>& electrolyte_ps,
        const std::string& filename, int save_interval);

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