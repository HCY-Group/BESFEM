// MaterialProperties.cpp
#include "../include/MaterialProperties.hpp"
#include "../include/Constants.hpp"
#include "mfem.hpp"
#include <cmath>
#include <fstream>

namespace MaterialProperties
{
    static double GetTableValues(double x, const mfem::Vector &ticks, const mfem::Vector &data)
    {
        const int n = ticks.Size();

        if (x <= ticks(0)){return data(0);}
        if (x >= ticks(n - 1)){return data(n - 1);}

        int idx = 0;
        while (idx < n - 2 && ticks(idx + 1) < x) { idx++;}

        const double dx = ticks(idx + 1) - ticks(idx);
        return data(idx) + (x - ticks(idx)) / dx * (data(idx + 1) - data(idx));
    }
    
    static double NMC_OCV(double c)
    {
        return (1.095 * c * c) - (8.234e-7 * std::exp(14.31 * c)) + (4.692 * std::exp(-0.5389 * c));
    }

    static double NMC_mu(double c)
    {
        double val = -Constants::Frd * NMC_OCV(c);
        return val;
    }

    static double NMC_i0(double c)
    {
        double val = -0.2 * (c - 0.37) - 1.559 - 0.9376 * std::tanh(8.961 * c - 3.195);
        return std::pow(10.0, val) * 1.0e-3;
    }

    static double NMC_diff(double c)
    {
        double val = (0.0277 - 0.084 * c + 0.1003 * c * c) * 1.0e-8;
        return val;
    }

    // static double LFP_i0(double c)
    // {
    //     double val = 1.2e-5 * std::pow(c, 0.35) * std::pow(1.0 - c, 2.5);
    //     return val * 6;
    // }

    static double LFP_i0(double c)
    {
        c = std::min(1.0 - 1.0e-8, std::max(1.0e-8, c));

        double i0_mA = 2.75
                    * (1.0 - std::exp(-18.0 * c))
                    * (1.0 - std::exp(-18.0 * (1.0 - c)));

        // std::cout << "LFP_i0 at c = " << c << " is " << i0_mA * 1.0e-6 << std::endl;

        return 4* i0_mA * 1.0e-6; // mA/cm^2 to A/cm^2

    }

    
    static double LFP_OCV(double c)
    {
        static mfem::Vector Ticks(201);
        static mfem::Vector chmPot(201);
        static bool loaded = false;

        if (!loaded)
        {
            std::ifstream myXfile("../inputs/materials/LFP_Chm_Pot_Ticks.txt");
            std::ifstream mydFfile("../inputs/materials/LFP_Chm_Pot.txt");

            if (!myXfile || !mydFfile)
            {
                mfem::mfem_error("Could not open LFP chemical potential input files.");
            }

            for (int i = 0; i < 201; i++) myXfile >> Ticks(i);
            for (int i = 0; i < 201; i++) mydFfile >> chmPot(i);

            loaded = true;
        }
        
        return (-1*GetTableValues(c, Ticks, chmPot)) + 3.4;    
    }
    
    static double LFP_mu(double c)
    {
	    double val = -Constants::Frd*LFP_OCV(c);
        // std::cout << "LFP_mu at c = " << c << " is " << val << std::endl;
        return val;
    }


    static double LFP_dmu_dc(double c)
    {
        const double h = 0.01;

        double c1 = std::max(0.0, c - h);
        double c2 = std::min(1.0, c + h);

        return (LFP_mu(c2) - LFP_mu(c1)) / (c2 - c1);
    }

    static double LFP_diff(double c)
    {
        // return 5.0e-14; // LFP diffusivity e-14 gives e-12 for mobility
        return 5.0e-12; // e-12 gives e-10 for mobility
    }

    // static double LFP_Mob(double c)
    // {
    //     c = std::min(1.0 - 1.0e-8, std::max(1.0e-8, c));

    //     double D = LFP_diff(c);
    //     double dmu_dc = LFP_dmu_dc(c);

    //     double M = D / std::abs(dmu_dc);

    //     if (!std::isfinite(M))
    //     {
    //         std::cout << "Bad mobility at c = " << c
    //                 << " dmu_dc = " << dmu_dc
    //                 << std::endl;
    //     }

    //     std::cout << "LFP_Mob at c = " << c << " is " << M << std::endl;

    //     return D / std::abs(dmu_dc);
    // }

    static double LFP_Mob(double c)
    {
        // c = std::min(1.0 - 1.0e-8, std::max(1.0e-8, c));

        // const double Mmin = 8.5e-13;
        // const double B    = 5.0e-14;

        // const double A1 = 4.5e-12;
        // const double A2 = 4.0e-12;

        // const double c1 = 0.07;
        // const double c2 = 0.86;

        // const double w1 = 0.04;
        // const double w2 = 0.04;

        // const double x1 = (c - c1) / w1;
        // const double x2 = (c - c2) / w2;

        // return Mmin
        //     + B * (c - 0.5) * (c - 0.5)
        //     + A1 / (1.0 + x1 * x1)
        //     + A2 / (1.0 + x2 * x2);

        return 1e-13;
    }

    double LFP_ChpValue(double c)
    {
        static mfem::Vector Ticks(201);
        static mfem::Vector chmPot(201);
        static bool loaded = false;

        if (!loaded)
        {
            std::ifstream myXfile("../inputs/materials/LFP_Chm_Pot_Ticks.txt");
            std::ifstream mydFfile("../inputs/materials/LFP_Chm_Pot.txt");

            if (!myXfile || !mydFfile)
            {
                mfem::mfem_error("Could not open LFP chemical potential input files.");
            }

            for (int i = 0; i < 201; i++) myXfile >> Ticks(i);
            for (int i = 0; i < 201; i++) mydFfile >> chmPot(i);

            loaded = true;
        }

        return GetTableValues(c, Ticks, chmPot) ;
    }

    static double Graphite_OCV(double c)
    {
        static mfem::Vector Ticks(101);
        static mfem::Vector OCV(101);
        static bool loaded = false;

        if (!loaded)
        {
            std::ifstream myXfile("../inputs/materials/C_Li_X_101.txt");
            std::ifstream myOCVfile("../inputs/materials/C_Li_O3_101.txt");

            if (!myXfile || !myOCVfile)
            {
                mfem::mfem_error("Could not open graphite OCV input files.");
            }

            for (int i = 0; i < 101; i++) myXfile >> Ticks(i);
            for (int i = 0; i < 101; i++) myOCVfile >> OCV(i);

            loaded = true;
        }

        return GetTableValues(c, Ticks, OCV);
    }

    double OCV(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::NMC:
                return NMC_OCV(c);

            case sim::MaterialType::LFP:
                return LFP_OCV(c);

            case sim::MaterialType::Graphite:
                return Graphite_OCV(c);

            default:
                mfem::mfem_error("Unknown material in OCV.");
                return 0.0;
        }
    }

    static double Graphite_i0(double c)
    {
        static mfem::Vector Ticks(101);
        static mfem::Vector i0(101);
        static bool loaded = false;

        if (!loaded)
        {
            std::ifstream myXfile("../inputs/materials/C_Li_X_101.txt");
            std::ifstream myi0file("../inputs/materials/C_Li_J2_101.txt");

            if (!myXfile || !myi0file)
            {
                mfem::mfem_error("Could not open graphite i0 input files.");
            }

            for (int i = 0; i < 101; i++) myXfile >> Ticks(i);
            for (int i = 0; i < 101; i++) myi0file >> i0(i);

            loaded = true;
        }

        return GetTableValues(c, Ticks, i0) * 1.0e-3; // mA/cm^2 to A/cm^2

    }

    double ExchangeCurrentDensity(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::NMC:
                return NMC_i0(c);

            case sim::MaterialType::LFP:
                return LFP_i0(c);

            case sim::MaterialType::Graphite:
                return Graphite_i0(c);

            default:
                mfem::mfem_error("Unknown material in ExchangeCurrentDensity.");
                return 0.0;
        }
    }

    static double Electrolyte_diff(double c)
    {
        return Constants::D0 * std::exp(-7.02 - 830 * c + 50000 * c * c);
    }

    double Diffusivity(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::NMC:
                return NMC_diff(c);

            case sim::MaterialType::Electrolyte:
                return Electrolyte_diff(c);

            case sim::MaterialType::LFP:
                return LFP_diff(c); // placeholder, constant diffusivity for LFP

            default:
                mfem::mfem_error("Material does not have a defined diffusivity.");
                return 0.0;
        }
    }

    static double Graphite_Mob(double c)
    {
        static mfem::Vector Ticks(101);
        static mfem::Vector Mob(101);
        static bool loaded = false;

        if (!loaded)
        {
            std::ifstream myXfile("../inputs/materials/C_Li_X_101.txt");
            std::ifstream mydFfile("../inputs/materials/C_Li_Mb5_101.txt");

            if (!myXfile || !mydFfile)
            {
                mfem::mfem_error("Could not open graphite mobility input files.");
            }

            for (int i = 0; i < 101; i++) myXfile >> Ticks(i);
            for (int i = 0; i < 101; i++) mydFfile >> Mob(i);
            for (int i = 0; i < 101; i++) Mob(i) *= 100.0 * 2.0/3.0;


            loaded = true;
        }
        
        return GetTableValues(c, Ticks, Mob);  
    }

    double Mobility(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::Graphite:
                return Graphite_Mob(c);

            case sim::MaterialType::LFP:
                return LFP_Mob(c);

            default:
                mfem::mfem_error("Unknown material in Mobility.");
                return 0.0;
        }
    }

    static double Graphite_mu(double c)
    {
        static mfem::Vector Ticks(101);
        static mfem::Vector chmPot(101);
        static bool loaded = false;

        if (!loaded)
        {
            std::ifstream myXfile("../inputs/materials/C_Li_X_101.txt");
            std::ifstream mydFfile("../inputs/materials/C_Li_M6_101.txt");

            if (!myXfile || !mydFfile)
            {
                mfem::mfem_error("Could not open graphite chemical potential input files.");
            }

            for (int i = 0; i < 101; i++) myXfile >> Ticks(i);
            for (int i = 0; i < 101; i++) mydFfile >> chmPot(i);

            loaded = true;
        }

        return GetTableValues(c, Ticks, chmPot);
    }

    double ChemicalPotential(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::Graphite:
                return Graphite_mu(c);
            
            case sim::MaterialType::NMC:
                return NMC_mu(c);

            case sim::MaterialType::LFP:
                return LFP_mu(c);
            default:
                mfem::mfem_error("Unknown material in Chemical Potential.");
                return 0.0;
        }
    }

    static double GraphiteConductivity(double c)
    {
        return 3.3; 
    }

    static double NMCConductivity(double c)
    {
        return (0.01929 + 0.7045 * tanh(2.399 * c) - 0.7238 * tanh(2.412 * c) - 4.2106e-6); 
    }

    static double LFPConductivity(double c)
    {
        // return 1e-11; // S/cm
        return 1e-4;
    }

    double Conductivity(sim::MaterialType material, double c)
    {
        switch (material)
        {
            case sim::MaterialType::Graphite:
                return GraphiteConductivity(c);

            case sim::MaterialType::NMC:
                return NMCConductivity(c);

            case sim::MaterialType::LFP:
                return LFPConductivity(c);

            default:
                mfem::mfem_error("Unknown material in Conductivity.");
                return 0.0;
        }
    }

    double SiteDensity(sim::MaterialType material)
    {
        switch (material)
        {
            case sim::MaterialType::Graphite:
                // std::cout << "using graphite density" << std::endl;
                return 0.0312;

            case sim::MaterialType::NMC:
                // std::cout << "using NMC density" << std::endl;
                return 0.0501;

            case sim::MaterialType::LFP:
                // std::cout << "using LFP density" << std::endl;
                return 0.02273544498;

            default:
                mfem::mfem_error("Unknown material in SiteDensity.");
                return 0.0;
        }
    }

}
