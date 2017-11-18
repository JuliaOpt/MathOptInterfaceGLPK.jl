using GLPK, Base.Test, MathOptInterface, MathOptInterfaceTests, MathOptInterfaceGLPK
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contlinear.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "intlinear.jl"))
include(joinpath(Pkg.dir("MathOptInterface"), "test", "contconic.jl"))

const MOIT = MathOptInterfaceTests
const MOIGLP = MathOptInterfaceGLPK

@testset "MathOptInterfaceGLPK" begin
    @testset "Linear tests" begin
        linconfig_nocertificate = MOIT.TestConfig(1e-8,1e-8,true,true,false)
        solver = MOIGLP.GLPKSolverLP()
        MOIT.contlineartest(solver, linconfig_nocertificate)
        
        linconfig_nocertificate_noduals = MOIT.TestConfig(1e-8,1e-8,true,false,false)
        solver_mip = MOIGLP.GLPKSolverMIP()
        MOIT.contlineartest(solver_mip, linconfig_nocertificate_noduals, ["linear8b","linear8c"])
    end

    @testset "Linear Conic tests" begin
        linconfig_nocertificate = MOIT.TestConfig(1e-8,1e-8,true,true,false)
        solver = MOIGLP.GLPKSolverLP()
        MOIT.lintest(solver, linconfig_nocertificate)
        
        linconfig_nocertificate_noduals = MOIT.TestConfig(1e-8,1e-8,true,false,false)
        solver_mip = MOIGLP.GLPKSolverMIP()
        MOIT.lintest(solver_mip, linconfig_nocertificate_noduals)
    end

    @testset "Integer Linear tests" begin
        intconfig = MOIT.TestConfig(1e-8,1e-8,true,true,true)
        solver = MOIGLP.GLPKSolverMIP()
        MOIT.intlineartest(solver, intconfig, ["int2","int1"])
    end
end
;