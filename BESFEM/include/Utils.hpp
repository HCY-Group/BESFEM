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

class Utils
{
public:

    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------
    Utils(Initialize_Geometry &geo, Domain_Parameters &para);

    // ------------------------------------------------------------------------
    // High-level functions (previously from Concentrations_Base)
    // ------------------------------------------------------------------------
    void SetInitialValue(mfem::ParGridFunction &Cn, double initial_value);
    void InitializeReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, double value);

    void InitializeReaction(mfem::ParGridFunction &Rx1, mfem::ParGridFunction &Rx2, mfem::ParGridFunction &Rx3, double value);

    void CalculateLithiation(mfem::ParGridFunction &Cn, mfem::ParGridFunction &psx, double gtps);

    void CalculateReactionInfx(mfem::ParGridFunction &Rx, double &xCrnt);

    void CalculateGlobalError(mfem::ParGridFunction &px0, mfem::ParGridFunction &potential, mfem::ParGridFunction &psx, double &globalerror, double gtPsx);

    double GetLithiation() const { return Xfr_; } ///< Return current degree of lithiation.

    // ------------------------------------------------------------------------
    // Directory builder
    // ------------------------------------------------------------------------
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

    // ------------------------------------------------------------------------
    // SaveSimulationSnapshot overloads
    // ------------------------------------------------------------------------
    static void SaveSimulationSnapshot(
        int t, const std::string &outdir,
        Initialize_Geometry &geometry,
        Domain_Parameters &domain_parameters,
        mfem::ParGridFunction &phA,
        mfem::ParGridFunction &phC,
        mfem::ParGridFunction &phE,
        mfem::ParGridFunction &CnA,
        mfem::ParGridFunction &CnC,
        mfem::ParGridFunction &CnE,
        mfem::ParGridFunction &CnApsi,
        mfem::ParGridFunction &CnCpsi,
        mfem::ParGridFunction &CnEpsi,
        mfem::ParGridFunction &CnP,
        int save_interval = 500);

    static void SaveSimulationSnapshot(
        int t, const std::string &outdir,
        Initialize_Geometry &geometry,
        Domain_Parameters &domain_parameters,
        mfem::ParGridFunction &phC,
        mfem::ParGridFunction &phE,
        mfem::ParGridFunction &CnC,
        mfem::ParGridFunction &CnE,
        mfem::ParGridFunction &CnCpsi,
        mfem::ParGridFunction &CnEpsi,
        int save_interval = 500);

private:

    // Saved geometry & domain info
    Initialize_Geometry &geometry_;
    Domain_Parameters   &domain_;

    // Cached quantities
    mfem::ParMesh *pmesh_;
    std::shared_ptr<mfem::ParFiniteElementSpace> fes_;

    int nE_, nC_, nV_;
    mfem::Vector EVol_;

    mfem::Vector EAvg_;
    mfem::Array<double> VtxVal_;

    mfem::ParGridFunction TmpF_;

    double Xfr_   = 0.0;   // lithiation
    double geCrnt_ = 0.0;  // global reaction
    double infx_   = 0.0;  // boundary flux
};
