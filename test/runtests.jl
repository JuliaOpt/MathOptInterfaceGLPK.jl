using GLPK, Base.Test, MathOptInterface, MathOptInterfaceGLPK
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contlinear.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "intlinear.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contconic.jl"))


# contlinear
linear1test(MathOptInterfaceGLPK.GLPKSolverLP())
linear2test(MathOptInterfaceGLPK.GLPKSolverLP())
linear3test(MathOptInterfaceGLPK.GLPKSolverLP())
linear4test(MathOptInterfaceGLPK.GLPKSolverLP())
linear5test(MathOptInterfaceGLPK.GLPKSolverLP())
linear6test(MathOptInterfaceGLPK.GLPKSolverLP())
linear7test(MathOptInterfaceGLPK.GLPKSolverLP())
linear10test(MathOptInterfaceGLPK.GLPKSolverLP())
linear8test(MathOptInterfaceGLPK.GLPKSolverLP()) # infeasible/unbounded
linear9test(MathOptInterfaceGLPK.GLPKSolverLP())
linear11test(MathOptInterfaceGLPK.GLPKSolverLP())

# # intlinear
# knapsacktest(MathOptInterfaceGLPK.GLPKSolverMIP())
# int3test(MathOptInterfaceGLPK.GLPKSolverMIP())
# int1test(MathOptInterfaceGLPK.GLPKSolverMIP())
# # int2test(MathOptInterfaceGLPK.GLPKSolverMIP()) # SOS

# # contconic
lin1tests(MathOptInterfaceGLPK.GLPKSolverLP())
lin2tests(MathOptInterfaceGLPK.GLPKSolverLP())
# lin3test(MathOptInterfaceGLPK.GLPKSolverLP()) # infeasible
# # lin4test(MathOptInterfaceGLPK.GLPKSolverLP()) # infeasible
