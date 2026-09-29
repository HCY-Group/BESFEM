/**
 * @file BESFEM_All.hpp
 * @brief Convenience header collecting the BESFEM simulation interfaces.
 */
#pragma once
#include "mfem.hpp"
#include "mpi.h"

#include "../include/Adjust.hpp"
#include "../include/BoundaryConditions.hpp"
#include "../include/Concentrations_Base.hpp"
#include "../include/Constants.hpp"
#include "../include/dist_solver.hpp"
#include "../include/Domain_Parameters.hpp"
#include "../include/ElectrodeCahnHilliard.hpp"
#include "../include/ElectrodeDiffusion.hpp"
#include "../include/ElectrodePotential.hpp"
#include "../include/ElectrolyteDiffusion.hpp"
#include "../include/ElectrolytePotential.hpp"
#include "../include/FEMOperators.hpp"
#include "../include/Initialize_Geometry.hpp"
#include "../include/MaterialProperties.hpp"
#include "../include/Potentials_Base.hpp"
#include "../include/Reaction.hpp"
#include "../include/readtiff.h"
#include "../include/SimTypes.hpp"
#include "../include/SimulationConfig.hpp"
#include "../include/SimulationState.hpp"
#include "../include/Utils.hpp"

