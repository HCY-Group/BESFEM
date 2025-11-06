#pragma once
#include <string>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <ctime>
#include "mfem.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/Domain_Parameters.hpp"

namespace Utils {

// ============================================================================
// Build timestamped output directory
// ============================================================================
inline std::string BuildRunOutdir(const char* mesh_file, int num_steps)
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

// ============================================================================
// Save all grid functions and mesh fields at a given timestep
// ============================================================================
inline void SaveSimulationSnapshot(
    int t,
    const std::string &outdir,
    Initialize_Geometry &geometry,
    Domain_Parameters &domain_parameters,
    mfem::ParGridFunction &phA_gf,
    mfem::ParGridFunction &phC_gf,
    mfem::ParGridFunction &phE_gf,
    mfem::ParGridFunction &CnA_gf,
    mfem::ParGridFunction &CnC_gf,
    mfem::ParGridFunction &CnE_gf,
    mfem::ParGridFunction &CnA_gf_psi,
    mfem::ParGridFunction &CnC_gf_psi,
    mfem::ParGridFunction &CnE_gf_psi,
    mfem::ParGridFunction &CnP_together,
    int save_interval = 500)
{
    if (t % save_interval != 0) return;

    std::ostringstream step_str;
    step_str << "_" << std::setw(5) << std::setfill('0') << t;
    std::string suffix = step_str.str();

    // Mesh and ψ fields
    geometry.parallelMesh->SaveAsOne((outdir + "/pmesh" + suffix).c_str());
    domain_parameters.pse->SaveAsOne((outdir + "/pse" + suffix).c_str());

    // Potentials
    phA_gf.SaveAsOne((outdir + "/phA" + suffix).c_str());
    phC_gf.SaveAsOne((outdir + "/phC" + suffix).c_str());
    phE_gf.SaveAsOne((outdir + "/phE" + suffix).c_str());

    // Raw concentrations
    CnA_gf.SaveAsOne((outdir + "/CnA_raw" + suffix).c_str());
    CnC_gf.SaveAsOne((outdir + "/CnC_raw" + suffix).c_str());
    CnE_gf.SaveAsOne((outdir + "/CnE_raw" + suffix).c_str());

    // ψ-weighted concentrations
    CnA_gf_psi = CnA_gf; CnA_gf_psi *= *domain_parameters.psA;
    CnA_gf_psi.SaveAsOne((outdir + "/CnA" + suffix).c_str());

    CnC_gf_psi = CnC_gf; CnC_gf_psi *= *domain_parameters.psC;
    CnC_gf_psi.SaveAsOne((outdir + "/CnC" + suffix).c_str());

    CnE_gf_psi = CnE_gf; CnE_gf_psi *= *domain_parameters.pse;
    CnE_gf_psi.SaveAsOne((outdir + "/CnE" + suffix).c_str());

    // Composite ψ-weighted concentration
    CnP_together = CnA_gf_psi;
    CnP_together += CnC_gf_psi;
    CnP_together.SaveAsOne((outdir + "/CnP" + suffix).c_str());
}

} // namespace besfem
