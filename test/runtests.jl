using GLPK, Base.Test, MathOptInterface, MathOptInterfaceTests, MathOptInterfaceGLPK

const MOIT = MathOptInterfaceTests
const MOIGLP = MathOptInterfaceGLPK

@testset "MathOptInterfaceGLPK" begin
    @testset "Linear tests" begin
        linconfig_nocertificate = MOIT.TestConfig(infeas_certificates=false)
        solverf() = GLPKSolverLPInstance()
        MOIT.contlineartest(solverf, linconfig_nocertificate)
        
        linconfig_nocertificate_noduals = MOIT.TestConfig(duals=false,infeas_certificates=false)
        solverf_mip() = GLPKSolverMIPInstance()
        MOIT.contlineartest(solverf_mip, linconfig_nocertificate_noduals, ["linear8b","linear8c"])
    end

    @testset "Linear Conic tests" begin
        linconfig_nocertificate = MOIT.TestConfig(infeas_certificates=false)
        solverf() = GLPKSolverLPInstance()
        MOIT.lintest(solverf, linconfig_nocertificate)
        
        linconfig_nocertificate_noduals = MOIT.TestConfig(duals=false,infeas_certificates=false)
        solverf_mip() = GLPKSolverMIPInstance()
        MOIT.lintest(solverf_mip, linconfig_nocertificate_noduals)
    end

    @testset "Integer Linear tests" begin
        intconfig = MOIT.TestConfig()
        solverf() = GLPKSolverMIPInstance()
        MOIT.intlineartest(solverf, intconfig, ["int2","int1"])
    end
end
;