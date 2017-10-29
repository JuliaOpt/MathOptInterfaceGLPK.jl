module MathOptInterfaceGLPK

import Base.show, Base.copy

# Standard LP interface
# importall MathProgBase.SolverInterface

using GLPK
using MathOptInterface
const MOI = MathOptInterface
using LinQuadOptInterface
const LQOI = LinQuadOptInterface


# export GLPKSolver
abstract type GLPKSolver <: LQOI.LinQuadSolver end
# export GLPKMIPSolver
# export GLPKLPSolver
struct GLPKSolverLP <: GLPKSolver
    presolve::Bool
    method::Symbol
    options
end
function GLPKSolverLP(presolve = false, method = :Simplex;kwargs...)
    method in [:Simplex, :Exact, :InteriorPoint] ||
    error("""
          Unknown method for GLPK LP solver: $method
                 Allowed methods:
                   :Simplex
                   :Exact
                   :InteriorPoint""")

    return GLPKSolverLP(presolve, method, kwargs)
end
struct GLPKSolverMIP <: GLPKSolver
    presolve::Bool
    options
end
function GLPKSolverMIP(presolve = false;kwargs...)
    return GLPKSolverMIP(presolve, kwargs)
end

import GLPK.Prob

const Model = GLPK.Prob
Model(env) = Model()

abstract type GLPKSolverInstance <: LQOI.LinQuadSolverInstance

mutable struct GLPKSolverLPInstance <: LGLPKSolverInstance

    LQOI.@LinQuadSolverInstanceBase

    method::Symbol
    param::Union{GLPK.SimplexParam, GLPK.InteriorParam}
   
end

function MOI.SolverInstance(s::GLPKSolverLP)
    if s.method == :Simplex || s.method == :Exact
        param = GLPK.SimplexParam()
        if s.presolve
            param.presolve = GLPK.ON
        end
    elseif s.method == :InteriorPoint
        param = GLPK.InteriorParam()
        if s.presolve
            warn("Ignored option: presolve")
        end
    else
        error("This is a bug")
    end
    param.msg_lev = GLPK.MSG_ERR
    for (k,v) in s.opts
        i = findfirst(x->x==k, fieldnames(typeof(param)))
        if i > 0
            t = typeof(param).types[i]
            setfield!(param, i, convert(t, v))
        else
            warn("Ignored option: $(string(k))")
        end
    end

    env = nothing
    m = GLPKSolverLPInstance(
        (LQOI.@LinQuadSolverInstanceBaseInit)...,
        s.method,
        param,
    )

    return m
end

mutable struct GLPKSolverMIPInstance <: LGLPKSolverInstance

    LQOI.@LinQuadSolverInstanceBase

    param::GLPK.IntoptParam
    smplxparam::GLPK.SimplexParam
    # lazycb::Union{Function,Void}
    # cutcb::Union{Function,Void}
    # heuristiccb::Union{Function,Void}
    # infocb::Union{Function,Void}
    objbound::Float64
    # cbdata::MathProgCallbackData
    binaries::Vector{Int}
    userlimit::Bool

    
end


function MOI.SolverInstance(s::GLPKSolverMIP)
    if s.method == :Simplex || s.method == :Exact
        param = GLPK.SimplexParam()
        if s.presolve
            param.presolve = GLPK.ON
        end
    elseif s.method == :InteriorPoint
        param = GLPK.InteriorParam()
        if s.presolve
            warn("Ignored option: presolve")
        end
    else
        error("This is a bug")
    end
    param.msg_lev = GLPK.MSG_ERR
    for (k,v) in s.opts
        i = findfirst(x->x==k, fieldnames(typeof(param)))
        if i > 0
            t = typeof(param).types[i]
            setfield!(param, i, convert(t, v))
        else
            warn("Ignored option: $(string(k))")
        end
    end

    env = nothing
    lpm = GLPKSolverLPInstance(
        (LQOI.@LinQuadSolverInstanceBaseInit)...,
        GLPK.IntoptParam(),
        GLPK.SimplexParam(),
        -Inf,
        Int[],
        false,
    )
    # lpm.cbdata = GLPKCallbackData(lpm)


    lpm.param.msg_lev = GLPK.MSG_ERR
    lpm.smplxparam.msg_lev = GLPK.MSG_ERR
    if s.presolve
        lpm.param.presolve = GLPK.ON
    end

    lpm.param.cb_func = cfunction(_internal_callback, Void, (Ptr{Void}, Ptr{Void}))
    lpm.param.cb_info = pointer_from_objref(lpm.cbdata)

    for (k,v) in s.opts
        if k in [:cb_func, :cb_info]
            warn("ignored option: $(string(k)); use the MathProgBase callback interface instead")
            continue
        end
        i = findfirst(x->x==k, fieldnames(typeof(lpm.param)))
        s = findfirst(x->x==k, fieldnames(typeof(lpm.smplxparam)))
        if !(i > 0 || s > 0)
            warn("Ignored option: $(string(k))")
            continue
        end
        if i > 0
            t = typeof(lpm.param).types[i]
            setfield!(lpm.param, i, convert(t, v))
        end
        if s > 0
            t = typeof(lpm.smplxparam).types[s]
            setfield!(lpm.smplxparam, s, convert(t, v))
        end
    end

    return lpm
end

#=
    inner wrapper
=#

#=
    Main
=#

# LinQuadSolver # Abstract type
# done above

# LQOI.lqs_setparam!(env, name, val)
# TODO fix this one
# LQOI.lqs_setparam!(m::GLPKSolverInstance, name, val) = GLPK.setfield!(m.inner, string(name), val)

# LQOI.lqs_setlogfile!(env, path)
# TODO fix this one
# LQOI.lqs_setlogfile!(m::GLPKSolverInstance, path) = GLPK.setlogfile(m.inner, path::String)

# LQOI.lqs_getprobtype(m)
# TODO - consider removing, apparently useless

#=
    Constraints
=#

cintvec(v::Vector) = convert(Vector{Int32}, v)
cdoublevec(v::Vector) = convert(Vector{Float64}, v)

# LQOI.lqs_chgbds!(m, colvec, valvec, sensevec)
LQOI.lqs_chgbds!(m::GLPK.Model, colvec, valvec, sensevec) = GLPK.cpx_chgbds!(m, colvec, valvec, sensevec)

# LQOI.lqs_getlb(m, col)
LQOI.lqs_getlb(m::GLPK.Model, col) = GLPK.cpx_getlb(m, col)
# LQOI.lqs_getub(m, col)
LQOI.lqs_getub(m::GLPK.Model, col) = GLPK.cpx_getub(m, col)

# LQOI.lqs_getnumrows(m)
LQOI.lqs_getnumrows(m::GLPK.Model) = GLPK.cpx_getnumrows(m)

# LQOI.lqs_addrows!(m, rowvec, colvec, coefvec, sensevec, rhsvec)
LQOI.lqs_addrows!(m::GLPK.Model, rowvec, colvec, coefvec, sensevec, rhsvec) = GLPK.cpx_addrows!(m::GLPK.Model, rowvec, colvec, coefvec, sensevec, rhsvec)

# LQOI.lqs_getrhs(m, rowvec)
LQOI.lqs_getrhs(m::GLPK.Model, row) = GLPK.cpx_getrhs(m, row)  

# colvec, coef = LQOI.lqs_getrows(m, rowvec)
# TODO improve
function LQOI.lqs_getrows(m::GLPK.Model, idx)
    return GLPK.cpx_getrows(m, idx)
end

# LQOI.lqs_getcoef(m, row, col) #??
# TODO improve
LQOI.lqs_getcoef(m::GLPK.Model, row, col) = GLPK.cpx_getcoef(m, row, col)

# LQOI.lqs_chgcoef!(m, row, col, coef)
# TODO SPLIT THIS ONE
LQOI.lqs_chgcoef!(m::GLPK.Model, row, col, coef)  = GLPK.cpx_chgcoef!(m::GLPK.Model, row, col, coef)

# LQOI.lqs_delrows!(m, row, row)
LQOI.lqs_delrows!(m::GLPK.Model, rowbeg, rowend) = GLPK.cpx_delrows!(m::GLPK.Model, rowbeg, rowend)

# LQOI.lqs_chgctype!(m, colvec, typevec)
# TODO fix types
LQOI.lqs_chgctype!(m::GLPK.Model, colvec, typevec) = GLPK.cpx_chgctype!(m::GLPK.Model, colvec, typevec)

# LQOI.lqs_chgsense!(m, rowvec, sensevec)
# TODO fix types
LQOI.lqs_chgsense!(m::GLPK.Model, rowvec, sensevec) = GLPK.cpx_chgsense!(m::GLPK.Model, rowvec, sensevec)

const VAR_TYPE_MAP = Dict{Symbol,Cchar}(
    :CONTINUOUS => Cchar('C'),
    :INTEGER => Cchar('I'),
    :BINARY => Cchar('B')
)
LQOI.lqs_vartype_map(m::GLPKSolverInstance) = VAR_TYPE_MAP

# LQOI.lqs_addsos(m, colvec, valvec, typ)
LQOI.lqs_addsos!(m::GLPK.Model, colvec, valvec, typ) = GLPK.add_sos!(m::GLPK.Model, typ, colvec, valvec)
# LQOI.lqs_delsos(m, idx, idx)
LQOI.lqs_delsos!(m::GLPK.Model, idx1, idx2) = GLPK.cpx_delsos!(m::GLPK.Model, idx1, idx2)

const SOS_TYPE_MAP = Dict{Symbol,Symbol}(
    :SOS1 => :SOS1,#Cchar('1'),
    :SOS2 => :SOS2#Cchar('2')
)
LQOI.lqs_sertype_map(m::GLPKSolverInstance) = SOS_TYPE_MAP

# LQOI.lqs_getsos(m, idx)
# TODO improve getting processes
function LQOI.lqs_getsos(m::GLPK.Model, idx)
    indices, weights, types = GLPK.cpx_getsos(m::GLPK.Model, idx)

    # types2 = Array{Symbol}(length(types))
    # for i in eachindex(types)
    #     if types[i] == Cchar('1')
    #         types2[i] = :SOS1
    #     elseif types[i] == Cchar('2')
    #         types2[i] = :SOS2
    #     end
    # end

    return indices, weights, types == Cchar('1') ? :SOS1 : :SOS2
end
# LQOI.lqs_getnumqconstrs(m)
LQOI.lqs_getnumqconstrs(m::GLPK.Model) = GLPK.cpx_getnumqconstrs(m::GLPK.Model)

# LQOI.lqs_addqconstr(m, cols,coefs,rhs,sense, I,J,V)
LQOI.lqs_addqconstr!(m::GLPK.Model, cols,coefs,rhs,sense, I,J,V) = GLPK.cpx_addqconstr!(m::GLPK.Model, cols,coefs,rhs,sense, I,J,V)

# LQOI.lqs_chgrngval
LQOI.lqs_chgrngval!(m::GLPK.Model, rows, vals) = GLPK.cpx_chgrngval!(m::GLPK.Model, rows, vals)

const CTR_TYPE_MAP = Dict{Symbol,Cchar}(
    :RANGE => Cchar('R'),
    :LOWER => Cchar('L'),
    :UPPER => Cchar('U'),
    :EQUALITY => Cchar('E')
)
LQOI.lqs_ctrtype_map(m::GLPKSolverInstance) = CTR_TYPE_MAP

#=
    Objective
=#

# LQOI.lqs_copyquad(m, intvec,intvec, floatvec) #?
LQOI.lqs_copyquad!(m::GLPK.Model, I, J, V) = GLPK.cpx_copyquad!(m::GLPK.Model, I, J, V)

# LQOI.lqs_chgobj(m, colvec,coefvec)
LQOI.lqs_chgobj!(m::GLPK.Model, colvec, coefvec)  = GLPK.cpx_chgobj!(m::GLPK.Model, colvec, coefvec) 

# LQOI.lqs_chgobjsen(m, symbol)
# TODO improve min max names
LQOI.lqs_chgobjsen!(m::GLPK.Model, symbol) = GLPK.cpx_chgobjsen!(m::GLPK.Model, symbol)
    

# LQOI.lqs_getobj(m)
LQOI.lqs_getobj(m::GLPK.Model) = GLPK.cpx_getobj(m::GLPK.Model) 

# lqs_getobjsen(m)
LQOI.lqs_getobjsen(m::GLPK.Model) = GLPK.cpx_getobjsen(m::GLPK.Model)

#=
    Variables
=#

# LQOI.lqs_getnumcols(m)
LQOI.lqs_getnumcols(m::GLPK.Model) = GLPK.cpx_getnumcols(m::GLPK.Model)

# LQOI.lqs_newcols!(m, int)
LQOI.lqs_newcols!(m::GLPK.Model, int) = GLPK.cpx_newcols!(m::GLPK.Model, int)

# LQOI.lqs_delcols!(m, col, col)
LQOI.lqs_delcols!(m::GLPK.Model, col, col2) = GLPK.cpx_delcols!(m::GLPK.Model, col, col2)

# LQOI.lqs_addmipstarts(m, colvec, valvec)
LQOI.lqs_addmipstarts!(m::GLPK.Model, colvec, valvec)  = GLPK.cpx_addmipstarts!(m::GLPK.Model, colvec, valvec) 

#=
    Solve
=#

# LQOI.lqs_mipopt!(m)
LQOI.lqs_mipopt!(m::GLPK.Model) = GLPK.cpx_mipopt!(m::GLPK.Model)

# LQOI.lqs_qpopt!(m)
LQOI.lqs_qpopt!(m::GLPK.Model) = GLPK.cpx_qpopt!(m::GLPK.Model)

# LQOI.lqs_lpopt!(m)
LQOI.lqs_lpopt!(m::GLPK.Model) = GLPK.cpx_lpopt!(m::GLPK.Model)


const TERMINATION_STATUS_MAP = Dict(
    GLPK.CPX_STAT_OPTIMAL                => MOI.Success,
    GLPK.CPX_STAT_UNBOUNDED              => MOI.UnboundedNoResult,
    GLPK.CPX_STAT_INFEASIBLE             => MOI.InfeasibleNoResult,
    GLPK.CPX_STAT_INForUNBD              => MOI.InfeasibleOrUnbounded,
    GLPK.CPX_STAT_OPTIMAL_INFEAS         => MOI.Success,
    GLPK.CPX_STAT_NUM_BEST               => MOI.NumericalError,
    GLPK.CPX_STAT_ABORT_IT_LIM           => MOI.IterationLimit,
    GLPK.CPX_STAT_ABORT_TIME_LIM         => MOI.TimeLimit,
    GLPK.CPX_STAT_ABORT_OBJ_LIM          => MOI.ObjectiveLimit,
    GLPK.CPX_STAT_ABORT_USER             => MOI.Interrupted,
    GLPK.CPX_STAT_OPTIMAL_FACE_UNBOUNDED => MOI.UnboundedNoResult,
    GLPK.CPX_STAT_ABORT_PRIM_OBJ_LIM     => MOI.ObjectiveLimit,
    GLPK.CPX_STAT_ABORT_DUAL_OBJ_LIM     => MOI.ObjectiveLimit,
    GLPK.CPXMIP_OPTIMAL                  => MOI.Success,
    GLPK.CPXMIP_OPTIMAL_TOL              => MOI.Success,
    GLPK.CPXMIP_INFEASIBLE               => MOI.InfeasibleNoResult,
    GLPK.CPXMIP_SOL_LIM                  => MOI.SolutionLimit,
    GLPK.CPXMIP_NODE_LIM_FEAS            => MOI.NodeLimit,
    GLPK.CPXMIP_NODE_LIM_INFEAS          => MOI.NodeLimit,
    GLPK.CPXMIP_TIME_LIM_FEAS            => MOI.TimeLimit,
    GLPK.CPXMIP_TIME_LIM_INFEAS          => MOI.TimeLimit,
    GLPK.CPXMIP_FAIL_FEAS                => MOI.OtherError,
    GLPK.CPXMIP_FAIL_INFEAS              => MOI.OtherError,
    GLPK.CPXMIP_MEM_LIM_FEAS             => MOI.MemoryLimit,
    GLPK.CPXMIP_MEM_LIM_INFEAS           => MOI.MemoryLimit,
    GLPK.CPXMIP_ABORT_FEAS               => MOI.Interrupted,
    GLPK.CPXMIP_ABORT_INFEAS             => MOI.Interrupted,
    GLPK.CPXMIP_OPTIMAL_INFEAS           => MOI.Success,
    GLPK.CPXMIP_FAIL_FEAS_NO_TREE        => MOI.MemoryLimit,
    GLPK.CPXMIP_FAIL_INFEAS_NO_TREE      => MOI.MemoryLimit,
    GLPK.CPXMIP_UNBOUNDED                => MOI.UnboundedNoResult,
    GLPK.CPXMIP_INForUNBD                => MOI.InfeasibleOrUnbounded
)

# LQOI.lqs_terminationstatus(m)
function LQOI.lqs_terminationstatus(model::GLPKSolverInstance)
    m = model.inner 

    code = GLPK.cpx_getstat(m)
    mthd, soltype, prifeas, dualfeas = GLPK.cpx_solninfo(m)

    
    if haskey(TERMINATION_STATUS_MAP, code)
        out = TERMINATION_STATUS_MAP[code]
        
        if code == GLPK.CPX_STAT_UNBOUNDED && prifeas > 0
            out = MOI.Success
        elseif code == GLPK.CPX_STAT_INFEASIBLE && dualfeas > 0
            out = MOI.Success
        end
        return out
    else
        error("Status $(code) has not been mapped to a MOI termination status.")
    end
end

function LQOI.lqs_primalstatus(model::GLPKSolverInstance)
    m = model.inner

    code = GLPK.cpx_getstat(m)
    mthd, soltype, prifeas, dualfeas = GLPK.cpx_solninfo(m)

    out = MOI.UnknownResultStatus

    if soltype in [GLPK.CPX_NONBASIC_SOLN, GLPK.CPX_BASIC_SOLN, GLPK.CPX_PRIMAL_SOLN]
        if prifeas > 0
            out = MOI.FeasiblePoint
        else
            out = MOI.InfeasiblePoint
        end
    end
    if code == GLPK.CPX_STAT_UNBOUNDED #&& prifeas > 0
        out = MOI.InfeasibilityCertificate
    end
    return out
end
function LQOI.lqs_dualstatus(model::GLPKSolverInstance)
    m = model.inner    

    code = GLPK.cpx_getstat(m)
    mthd, soltype, prifeas, dualfeas = GLPK.cpx_solninfo(m)
    if !LQOI.hasinteger(model)
        if soltype in [GLPK.CPX_NONBASIC_SOLN, GLPK.CPX_BASIC_SOLN]
            if dualfeas > 0
                out = MOI.FeasiblePoint
            else
                out = MOI.InfeasiblePoint
            end
        else
            out = MOI.UnknownResultStatus
        end
        if code == GLPK.CPX_STAT_INFEASIBLE && dualfeas > 0
            out = MOI.InfeasibilityCertificate
        end
        return out
    end
    return MOI.UnknownResultStatus
end


# LQOI.lqs_getx!(m, place)
LQOI.lqs_getx!(m::GLPK.Model, place) = GLPK.cpx_getx!(m::GLPK.Model, place) 

# LQOI.lqs_getax!(m, place)
LQOI.lqs_getax!(m::GLPK.Model, place) = GLPK.cpx_getax!(m::GLPK.Model, place)

# LQOI.lqs_getdj!(m, place)
LQOI.lqs_getdj!(m::GLPK.Model, place) = GLPK.cpx_getdj!(m::GLPK.Model, place)

# LQOI.lqs_getpi!(m, place)
LQOI.lqs_getpi!(m::GLPK.Model, place) = GLPK.cpx_getpi!(m::GLPK.Model, place)

# LQOI.lqs_getobjval(m)
LQOI.lqs_getobjval(m::GLPK.Model) = GLPK.cpx_getobjval(m::GLPK.Model)

# LQOI.lqs_getbestobjval(m)
LQOI.lqs_getbestobjval(m::GLPK.Model) = GLPK.cpx_getbestobjval(m::GLPK.Model)

# LQOI.lqs_getmiprelgap(m)
LQOI.lqs_getmiprelgap(m::GLPK.Model) = GLPK.cpx_getmiprelgap(m::GLPK.Model)

# LQOI.lqs_getitcnt(m)
LQOI.lqs_getitcnt(m::GLPK.Model)  = GLPK.cpx_getitcnt(m::GLPK.Model)

# LQOI.lqs_getbaritcnt(m)
LQOI.lqs_getbaritcnt(m::GLPK.Model) = GLPK.cpx_getbaritcnt(m::GLPK.Model)

# LQOI.lqs_getnodecnt(m)
LQOI.lqs_getnodecnt(m::GLPK.Model) = GLPK.cpx_getnodecnt(m::GLPK.Model)

# LQOI.lqs_dualfarkas(m, place)
LQOI.lqs_dualfarkas!(m::GLPK.Model, place) = GLPK.cpx_dualfarkas!(m::GLPK.Model, place)

# LQOI.lqs_getray(m, place)
LQOI.lqs_getray!(m::GLPK.Model, place) = GLPK.cpx_getray!(m::GLPK.Model, place)


MOI.free!(m::GLPKSolverInstance) = GLPK.free_model(m.inner)

"""
    writeproblem(m::AbstractSolverInstance, filename::String)
Writes the current problem data to the given file.
Supported file types are solver-dependent.
"""
writeproblem(m::GLPKSolverInstance, filename::String, flags::String="") = GLPK.write_model(m.inner, filename)


LQOI.lqs_make_problem_type_continuous(m::GLPK.Model) = GLPK._make_problem_type_continuous(m)
end # module