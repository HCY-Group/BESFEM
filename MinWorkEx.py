import mfem.par as mfem

# test initializing
mesh = mfem.Mesh()
print('mesh initialized successfully')
pargf = mfem.ParGridFunction()
print('pargridfunction initialized successfully')

distsolver = mfem.dist_solver.HeatDistanceSolver(1.0)
print('distsolver initialized successfully')


