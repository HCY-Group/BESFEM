#include "mfem.hpp"

using namespace mfem;
using namespace std;

double BvP;

void InitializePotP(mfem::ParFiniteElementSpace &fespace)
{
    ParGridFunction phP(&fespace);  // Electropotential in particle
    ParGridFunction kap(&fespace);  // Conductivity in particle
    ParGridFunction RpP(&fespace);  // Reaction
    ParGridFunction pP0(&fespace);  // Initial potential

    BvP = 2.9395;
    phP = BvP;

    // Stiffness matrix
    HypreParMatrix KmP;
    ParBilinearForm *Kp2;

    // Solver setup
    CGSolver cgPP(MPI_COMM_WORLD);
    cgPP.SetRelTol(1e-7);
    cgPP.SetMaxIter(200);

    // Force vector
    ParLinearForm *Bp2;
    ParLinearForm Fpt(&fespace);
    HypreParVector Fpb(&fespace);

    // Initial solution vector
    HypreParVector Xs0(&fespace);

    //cout << "BvPTEST: " << BvP << std::endl;
}
