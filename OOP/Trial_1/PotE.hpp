#ifndef POTE_HPP
#define POTE_HPP

#include "mfem.hpp"
#include "PotP.hpp"
#include "MeshHandler.hpp"
#include "CnE.hpp"
#include "CnP.hpp"

#include <memory>

class CnE;

class PotE {

public: 

    PotE(MeshHandler &mesh_handler, PotP &potp, CnE &cne);

    void Initialize();
    void TimeStep(double dt);

    double Gettc1() const { return tc1; }
    double Gettc2() const { return tc2; }

    const mfem::ParGridFunction& GetphE() const { return phE; }
    double GetBvE() const { return BvE; } 




private:
    MeshHandler &mesh_handler;
    PotP &potp;
    CnE &cne;


    mfem::ParFiniteElementSpace *fespace;
    mfem::ParGridFunction CnEGridFunction;
    // std::unique_ptr<mfem::ParMesh> pmesh;
    mfem::ParMesh *pmesh;

    mfem::ParGridFunction phE; // electropot in electrolyte
    mfem::ParGridFunction kpl; // electrolyte conductivity
    mfem::ParGridFunction RpE; // reaction rate for electrolyte
    mfem::ParGridFunction pE0;
    mfem::ParGridFunction Dmp;

    mfem::GridFunctionCoefficient cDm;

    mfem::HypreParMatrix Kml;
    mfem::HypreParMatrix Kdm;

    mfem::ParLinearForm B1t;
    mfem::HypreParVector X1v, B1v, CeVn, LpCe;

    double BvP;
    double BvE;

    double tc1;
    double tc2;

    double Vcell;
    double dffe;

};

#endif // POTE_HPP
