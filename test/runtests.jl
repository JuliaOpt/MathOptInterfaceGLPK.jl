using GLPK, Base.Test, MathOptInterface, MathOptInterfaceGLPK
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contlinear.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "intlinear.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contconic.jl"))
# include(joinpath(Pkg.dir("MathOptInterface"), "test", "contquadratic.jl"))


# contlinear
linear1test(MathOptInterfaceGLPK.GLPKSolverMIP())
# linear2test(MathOptInterfaceGLPK.GLPKSolver())
# linear3test(MathOptInterfaceGLPK.GLPKSolver())
# linear4test(MathOptInterfaceGLPK.GLPKSolver())
# linear5test(MathOptInterfaceGLPK.GLPKSolver())
# linear6test(MathOptInterfaceGLPK.GLPKSolver())
# linear7test(MathOptInterfaceGLPK.GLPKSolver())
# # linear8test(MathOptInterfaceGLPK.GLPKSolver()) # infeasible/unbounded
# linear9test(MathOptInterfaceGLPK.GLPKSolver())
# linear10test(MathOptInterfaceGLPK.GLPKSolver())
# linear11test(MathOptInterfaceGLPK.GLPKSolver())

# # intlinear
# knapsacktest(MathOptInterfaceGLPK.GLPKSolver())
# int1test(MathOptInterfaceGLPK.GLPKSolver())
# # int2test(MathOptInterfaceGLPK.GLPKSolver()) # SOS
# int3test(MathOptInterfaceGLPK.GLPKSolver())

# # contconic
# lin1tests(MathOptInterfaceGLPK.GLPKSolver())
# lin2tests(MathOptInterfaceGLPK.GLPKSolver())
# # lin3test(MathOptInterfaceGLPK.GLPKSolver()) # infeasible
# # lin4test(MathOptInterfaceGLPK.GLPKSolver()) # infeasible
