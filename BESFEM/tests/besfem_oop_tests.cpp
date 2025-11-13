// tests/besfem_oop_tests.cpp

#include "mfem.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <cstdio>

#include "../inputs/Constants.hpp"
#include "../include/Initialize_Geometry.hpp"

using std::cout;
using std::endl;

// ----------------- check & info statments -----------------
static int fails = 0;

#define CHECK(cond, msg) \
    do { \
        if (!(cond)) { \
            std::cerr << "[FAIL] " << msg << std::endl; \
            ++fails; \
        } \
    } while (0)

#define INFO(msg) \
    do { \
        std::cout << "[INFO] " << msg << std::endl; \
    } while (0)

// ----------------- helpers -----------------

// returns true if path exists and is a regular file
static bool FileExists(const std::string &path)
{
    struct stat st;
    return ::stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

// Write a simple distance file: one double per line
static void WriteDistanceFile(const std::string &path,
                              const std::vector<double> &vals)
{
    std::ofstream ofs(path);
    if (!ofs.is_open()) {
        std::cerr << "[FAIL] Could not open " << path << " for writing\n";
        ++fails;
        return;
    }
    ofs.setf(std::ios::scientific);
    ofs.precision(16);
    for (double v : vals) {
        ofs << v << "\n";
    }
}

// Read doubles from a file
static std::vector<double> ReadDistanceFile(const std::string &path)
{
    std::ifstream ifs(path);
    std::vector<double> vals;
    if (!ifs.is_open()) {
        std::cerr << "[FAIL] Could not open " << path << " for reading\n";
        ++fails;
        return vals;
    }
    double x;
    while (ifs >> x) {
        vals.push_back(x);
    }
    return vals;
}

// ----------------- tests -----------------

// check if all |v| <= 1, then file should be unchanged and no backup created
static void Test_AdjustDistanceFile_NoScaling()
{
    INFO("Test_AdjustDistanceFile_NoScaling");

    const std::string fname  = "dist_noscale.dsF";
    const std::string backup = fname + ".orig";

    // Clean up any leftovers
    std::remove(fname.c_str());
    std::remove(backup.c_str());

    // Original values: all within [-1, 1]
    std::vector<double> original = {0.5, -0.8, 0.1};
    WriteDistanceFile(fname, original);

    Initialize_Geometry geom;
    geom.AdjustDistanceFile(fname.c_str());

    // File should be unchanged
    auto after = ReadDistanceFile(fname);
    CHECK(after.size() == original.size(), "NoScaling: expected size " + std::to_string(original.size()) + ", got " + std::to_string(after.size()));
    if (after.size() == original.size())
    {
        for (size_t i = 0; i < original.size(); ++i)
        {
            CHECK(after[i] == original[i], "NoScaling: value at index " + std::to_string(i) + " should be unchanged");
        }
    }

    // No backup file should exist
    CHECK(!FileExists(backup), "NoScaling: backup file should NOT exist");

    // Cleanup
    std::remove(fname.c_str());
    std::remove(backup.c_str());
}

// check if max |v| > 1, then backup is written and file scaled by Constants::dh
static void Test_AdjustDistanceFile_ScalingAndBackup()
{
    INFO("Test_AdjustDistanceFile_ScalingAndBackup");

    const std::string fname  = "dist_scale.dsF";
    const std::string backup = fname + ".orig";

    // Clean up any leftovers
    std::remove(fname.c_str());
    std::remove(backup.c_str());

    // Some values with |v| > 1
    std::vector<double> original = {-2.0, 0.0, 3.0};
    WriteDistanceFile(fname, original);

    Initialize_Geometry geom;
    geom.AdjustDistanceFile(fname.c_str());

    // Backup file must exist and contain original values
    CHECK(FileExists(backup), "Scaling: backup file should exist");
    if (FileExists(backup))
    {
        auto backup_vals = ReadDistanceFile(backup);
        CHECK(backup_vals.size() == original.size(), "Scaling: backup size matches original");
        if (backup_vals.size() == original.size())
        {
            for (size_t i = 0; i < original.size(); ++i)
            {
                CHECK(backup_vals[i] == original[i], "Scaling: backup value at index " + std::to_string(i) + " matches original");
            }
        }
    }

    // Main file should be scaled by Constants::dh
    auto scaled_vals = ReadDistanceFile(fname);
    CHECK(scaled_vals.size() == original.size(), "Scaling: scaled size matches original");
    if (scaled_vals.size() == original.size())
    {
        for (size_t i = 0; i < original.size(); ++i)
        {
            double expected = original[i] * Constants::dh;
            double diff = std::fabs(scaled_vals[i] - expected);
            CHECK(diff < 1e-12, "Scaling: index " + std::to_string(i) + " expected " + std::to_string(expected) + " got " + std::to_string(scaled_vals[i]));
        }
    }

    // Cleanup
    std::remove(fname.c_str());
    std::remove(backup.c_str());
}

// ----------------- main -----------------

int main(int argc, char **argv)
{
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();

    Test_AdjustDistanceFile_NoScaling();
    Test_AdjustDistanceFile_ScalingAndBackup();

    if (fails == 0) {
        cout << "\nALL TESTS PASSED \n";
    } else {
        cout << "\n" << fails << " TEST(S) FAILED \n";
    }

    mfem::Hypre::Finalize();
    mfem::Mpi::Finalize();
    return (fails == 0) ? 0 : 1;
}
