using GLPK, Base.Test, MathOptInterface, MathOptInterface.Test, MathOptInterfaceGLPK

const MOIT = MathOptInterface.Test
const MOIGLP = MathOptInterfaceGLPK

@testset "MathOptInterfaceGLPK" begin
    @testset "Linear tests" begin
        linconfig_nocertificate = MOIT.TestConfig(infeas_certificates=false)
        solver = GLPKOptimizerLP()
        MOIT.contlineartest(solver, linconfig_nocertificate)
        
        linconfig_nocertificate_noduals = MOIT.TestConfig(duals=false,infeas_certificates=false)
        solver_mip = GLPKOptimizerMIP()
        MOIT.contlineartest(solver_mip, linconfig_nocertificate_noduals, ["linear8b","linear8c"])
    end

    @testset "Linear Conic tests" begin
        linconfig_nocertificate = MOIT.TestConfig(infeas_certificates=false)
        solver = GLPKOptimizerLP()
        MOIT.lintest(solver, linconfig_nocertificate)
        
        linconfig_nocertificate_noduals = MOIT.TestConfig(duals=false,infeas_certificates=false)
        solver_mip = GLPKOptimizerMIP()
        MOIT.lintest(solver_mip, linconfig_nocertificate_noduals)
    end

    @testset "Integer Linear tests" begin
        intconfig = MOIT.TestConfig()
        solver_mip = GLPKOptimizerMIP()
        MOIT.intlineartest(solver_mip, intconfig, ["int2","int1"])
    end
end
;