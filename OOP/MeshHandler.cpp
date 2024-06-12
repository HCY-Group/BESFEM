#include "MeshHandler.hpp"
#include "mfem.hpp"
#include "mpi.h"
#include <cmath>
#include <memory>

using namespace mfem;
using namespace std;

MeshHandler::MeshHandler(const char* mesh_file, const char* dsF_file, int order)
    : mesh_file(mesh_file), dsF_file(dsF_file), order(order), dh(0.2e-4), zeta(0.375), eps(1.0e-6), rho(0.0501), Cr(3.0),
      gmesh(mesh_file), gtPsi(0.0), gtPse(0.0), gTrgI(0.0) {}

void MeshHandler::InitializeMesh()
{
    gmesh.EnsureNCMesh(true);
    gFespace = make_unique<FiniteElementSpace>(&gmesh, new H1_FECollection(order, gmesh.Dimension()));

    ReadGlobalDistanceFunction();

    pmesh = make_unique<ParMesh>(MPI_COMM_WORLD, gmesh);
    fespace = make_unique<ParFiniteElementSpace>(pmesh.get(), new H1_FECollection(order, pmesh->Dimension()));
    dsF = make_unique<ParGridFunction>(fespace.get());

    // Map local distance function from global one
    Array<HYPRE_BigInt> E_L2G;
    pmesh->GetGlobalElementIndices(E_L2G);

    Array<int> gVTX(pow(2, pmesh->Dimension()));
    Array<int> VTX(pow(2, pmesh->Dimension()));
    for (int ei = 0; ei < pmesh->GetNE(); ei++) {
        int gei = E_L2G[ei];
        gmesh.GetElementVertices(gei, gVTX);
        pmesh->GetElementVertices(ei, VTX);
        for (int vi = 0; vi < pow(2, pmesh->Dimension()); vi++) {
            (*dsF)(VTX[vi]) = (*gDsF)(gVTX[vi]);
        }
    }

    InterpolateDomainParameters();
    CalculateTotalPsi();
    CalculateTotalPse();
    //CalculateTargetCurrent();
}

void MeshHandler::ReadGlobalDistanceFunction()
{
    gDsF = make_unique<GridFunction>(gFespace.get());
    ifstream myfile(dsF_file);
    for (int gi = 0; gi < gDsF->Size(); gi++) {
        myfile >> (*gDsF)(gi);
    }
    myfile.close();
}

void MeshHandler::InterpolateDomainParameters()
{
    psi = make_unique<ParGridFunction>(fespace.get());
    pse = make_unique<ParGridFunction>(fespace.get());
    AvP = make_unique<ParGridFunction>(fespace.get());

    for (int vi = 0; vi < fespace->GetNV(); vi++) {
        (*psi)(vi) = 0.5 * (1.0 + tanh((*dsF)(vi) / (zeta * dh)));
        (*pse)(vi) = 1.0 - (*psi)(vi);
        (*AvP)(vi) = -(pow(tanh((*dsF)(vi) / (zeta * dh)), 2) - 1.0) / (2 * zeta * dh);

        if ((*psi)(vi) < eps) { (*psi)(vi) = eps; }
        if ((*pse)(vi) < eps) { (*pse)(vi) = eps; }
    }

    for (int vi = 0; vi < fespace->GetNV(); vi++) {
        if ((*AvP)(vi) * dh < 1.0e-3) { (*AvP)(vi) = 0.0; }
    }
}

void MeshHandler::CalculateTotalPsi()
{
    double tPsi = 0.0;
    Vector EVol(pmesh->GetNE());
    for (int ei = 0; ei < pmesh->GetNE(); ei++) {
        EVol(ei) = pmesh->GetElementVolume(ei);
    }

    Vector EAvg(pmesh->GetNE());
    for (int ei = 0; ei < pmesh->GetNE(); ei++) {
        Array<double> VtxVal(pow(2, pmesh->Dimension()));
        psi->GetNodalValues(ei, VtxVal);
        double val = 0.0;
        for (int vt = 0; vt < pow(2, pmesh->Dimension()); vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / pow(2, pmesh->Dimension());
        tPsi += EAvg(ei) * EVol(ei);
    }

    MPI_Allreduce(&tPsi, &gtPsi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    CalculateTargetCurrent(tPsi);
}

void MeshHandler::CalculateTotalPse()
{
    double tPse = 0.0;
    Vector EVol(pmesh->GetNE());
    for (int ei = 0; ei < pmesh->GetNE(); ei++) {
        EVol(ei) = pmesh->GetElementVolume(ei);
    }

    Vector EAvg(pmesh->GetNE());
    for (int ei = 0; ei < pmesh->GetNE(); ei++) {
        Array<double> VtxVal(pow(2, pmesh->Dimension()));
        pse->GetNodalValues(ei, VtxVal);
        double val = 0.0;
        for (int vt = 0; vt < pow(2, pmesh->Dimension()); vt++) {
            val += VtxVal[vt];
        }
        EAvg(ei) = val / pow(2, pmesh->Dimension());
        tPse += EAvg(ei) * EVol(ei);
    }

    MPI_Allreduce(&tPse, &gtPse, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MeshHandler::CalculateTargetCurrent(double tPsi)
{
    double trgI = tPsi * rho * (0.9 - 0.3) / (3600.0 / Cr);
    MPI_Allreduce(&trgI, &gTrgI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MeshHandler::PrintMeshInfo()
{
    cout << "Number of vertices: " << pmesh->GetNV() << endl;
    cout << "Number of elements: " << pmesh->GetNE() << endl;

    Vector Rmin, Rmax;
    gmesh.GetBoundingBox(Rmin, Rmax);
    double L_w = Rmax(1) - Rmin(1);
    cout << "West boundary size: " << L_w << endl;

    cout << "Total Psi: " << gtPsi << endl;
    cout << "Total Pse: " << gtPse << endl;
    cout << "Target Current: " << gTrgI << endl;
}
