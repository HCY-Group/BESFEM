/**
 * @file MaterialProperties.hpp
 * @brief Material-property functions used throughout BESFEM.
 *
 * This file provides material-dependent electrochemical properties, including
 * open-circuit voltage, chemical potential, diffusivity, mobility,
 * conductivity, exchange current density, and site density.
 */

#pragma once
#include "SimTypes.hpp"
#include <string>
#include <unordered_map>

/**
 * @namespace MaterialProperties
 * @brief Material-dependent electrochemical property models.
 *
 * The functions in this namespace provide constitutive relationships for each
 * supported electrode material. These properties are used by the reaction,
 * concentration, and potential solvers.
 */
namespace MaterialProperties
{
    /** Load per-material default text files, then apply material.* overrides.
     * Override paths and materials_dir are relative to the config file.
     * Replaces the active database atomically; throws std::runtime_error for
     * invalid input. Call before validation/initialization, not during a solve.
     * The existing material API supports one active configuration per process.
     */
    void Configure(const std::unordered_map<std::string, std::string>& values,
                   const std::string& config_file);

    /**
     * @brief Return the open-circuit voltage (OCV).
     *
     * @param material Electrode material.
     * @param c Normalized lithium concentration.
     * @return Open-circuit voltage (V).
     */
    double OCV(sim::MaterialType material, double c);

    /**
     * @brief Return the chemical potential.
     *
     * @param material Electrode material.
     * @param c Normalized lithium concentration.
     * @return Chemical potential.
     */
    double ChemicalPotential(sim::MaterialType material, double c);

    /**
     * @brief Return the exchange current density.
     *
     * @param material Electrode material.
     * @param c Normalized lithium concentration.
     * @return Exchange current density.
     */
    double ExchangeCurrentDensity(sim::MaterialType material, double c);

    /**
     * @brief Return the tabulated LFP chemical-potential value.
     *
     * This function evaluates the tabulated chemical-potential curve used for
     * lithium iron phosphate (LFP).
     *
     * @param c Normalized lithium concentration.
     * @return LFP chemical-potential value.
     */
    double LFP_ChpValue(double c);

    /**
     * @brief Return the lithium diffusivity.
     *
     * @param material Electrode material.
     * @param c Normalized lithium concentration.
     * @return Diffusivity.
     */
    double Diffusivity(sim::MaterialType material, double c);

    /**
     * @brief Return the mobility.
     *
     * @param material Electrode material.
     * @param c Normalized lithium concentration.
     * @return Mobility.
     */
    double Mobility(sim::MaterialType material, double c);

    /**
     * @brief Return the electrical conductivity.
     *
     * @param material Electrode material.
     * @param c Normalized lithium concentration.
     * @return Electrical conductivity.
     */
    double Conductivity(sim::MaterialType material, double c);

    /**
     * @brief Return the lithium site density.
     *
     * @param material Electrode material.
     * @return Site density.
     */
    double SiteDensity(sim::MaterialType material);
}