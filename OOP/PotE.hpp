#ifndef POTE_HPP
#define POTE_HPP

#include "mfem.hpp"
#include "PotP.hpp"
#include "MeshHandler.hpp"

#include <memory>

class PotE {

public: 

    PotE(MeshHandler &mesh_handler, PotP &potp);

    void Initialize();

    double Gettc1() const { return tc1; }
    double Gettc2() const { return tc2; }




private:
    MeshHandler &mesh_handler;
    PotP &potp;


    mfem::ParFiniteElementSpace *fespace;
    mfem::ParGridFunction phE; // electropot in electrolyte
    mfem::ParGridFunction kpl; // electrolyte conductivity
    mfem::ParGridFunction RpE; // reaction rate for electrolyte
    mfem::ParGridFunction pE0;
    mfem::ParGridFunction Dmp;

    mfem::HypreParMatrix Kml;
    mfem::HypreParMatrix Kdm;

    double BvP;
    double BvE;

    double tc1;
    double tc2;

    double Vcell;
    double dffe;

};

#endif // POTE_HPP
