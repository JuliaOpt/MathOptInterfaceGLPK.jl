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
function GLPKSolverMIP(presolve = false; kwargs...)
    return GLPKSolverMIP(presolve, kwargs)
end

import GLPK.Prob

const Model = GLPK.Prob
Model(env) = Model()

abstract type GLPKSolverInstance <: LQOI.LinQuadSolverInstance end

type GLPKSolverLPInstance <: GLPKSolverInstance

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

mutable struct GLPKSolverMIPInstance <: GLPKSolverInstance

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

    env = nothing
    lpm = GLPKSolverMIPInstance(
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

    # lpm.param.cb_func = cfunction(_internal_callback, Void, (Ptr{Void}, Ptr{Void}))
    # lpm.param.cb_info = pointer_from_objref(lpm.cbdata)

    for (k,v) in s.options
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
function LQOI.lqs_chgbds!(m::GLPK.Prob, colvec, valvec, sensevec) 
    colub = Inf
    collb = -Inf
    bt = GLPK.DB
    for i in eachindex(colvec)
        if sensevec[i] == Cchar('E')
            colub = valvec[i]
            collb = valvec[i]
            bt = GLPK.FX
        elseif sensevec[i] == Cchar('G')
            collb = valvec[i]
            colub = Inf
            u = GLPK.get_col_ub(m, colvec[i])
            if u < Inf
                bt = GLPK.DB
                colub = u
            else
                bt = GLPK.LO
            end
        elseif sensevec[i] == Cchar('L')
            colub = valvec[i]
            collb = -Inf
            l = GLPK.get_col_lb(m, colvec[i])
            if l > -Inf
                bt = GLPK.DB
                collb = l
            else
                bt = GLPK.UP
            end
        else
            error("invalid bound type")
        end
        if colub == Inf && collb == -Inf
            bt = GLPK.FR
        end
        GLPK.set_col_bnds(m, colvec[i], bt, collb, colub)
    end 

end


# LQOI.lqs_getlb(m, col)
LQOI.lqs_getlb(m::GLPK.Prob, col) = GLPK.get_col_lb(m, col)

# LQOI.lqs_getub(m, col)
LQOI.lqs_getub(m::GLPK.Prob, col) = GLPK.get_col_ub(m, col)

# LQOI.lqs_getnumrows(m)
LQOI.lqs_getnumrows(m::GLPK.Prob) = GLPK.get_num_rows(m)

# LQOI.lqs_addrows!(m, rowvec, colvec, coefvec, sensevec, rhsvec)
function LQOI.lqs_addrows!(m::GLPK.Prob, rowvec, colvec, coefvec, sensevec, rhsvec) 
    
    nrows = length(rhsvec)

    if nrows <= 0
        error("nor row to be added")
    elseif nrows == 1
        addrow!(m::GLPK.Prob, colvec::Vector, coefvec::Vector, sensevec[1], rhsvec[1])

    else

        colidx = [Int[] for i in 1:nrows]
        colval = [Float64[] for i in 1:nrows]

        for i in 1:nrows
            nels = count(x->x==i, rowvec)
            sizehint!(colidx, nels)
            sizehint!(colval, nels)
        end

        for i in eachindex(rowvec)
            push!(colidx[rowvec[i]], colvec[i])
            push!(colval[rowvec[i]], coefvec[i])
        end

        for i in 1:nrows
            addrow!(lp::GLPK.Prob, colidx[i], colval[i]::Vector, sense[i]::Cchar, rhs[i]::Real)
        end
        
    end
    nothing
end
function addrow!(lp::GLPK.Prob, colidx::Vector, colcoef::Vector, sense::Cchar, rhs::Real)
    if length(colidx) != length(colcoef)
        error("colidx and colcoef have different legths")
    end
    GLPK.add_rows(lp, 1)
    m = GLPK.get_num_rows(lp)
    GLPK.set_mat_row(lp, m, colidx, colcoef)
    
    if sense == Cchar('E')
        bt = GLPK.FX
        rowlb = rhs
        rowub = rhs
    elseif sense == Cchar('G')
        bt = GLPK.LO
        rowlb = rhs
        rowub = Inf
    elseif sense == Cchar('L') 
        bt = GLPK.UP
        rowlb = -Inf
        rowub = rhs
    else
        error("row type not valid")
        bt = GLPK.FR
    end
    GLPK.set_row_bnds(lp, m, bt, rowlb, rowub)
    return
end


# LQOI.lqs_getrhs(m, rowvec)
function LQOI.lqs_getrhs(m::GLPK.Prob, row)
    sense = GLPK.get_row_type(m, row)
    if sense == GLPK.LO
        return GLPK.get_row_lb(m, row)
    elseif sense == GLPK.FX
        return GLPK.get_row_lb(m, row)
    elseif sense == GLPK.DB
        return GLPK.get_row_lb(m, row)
    else
        return GLPK.get_row_ub(m, row)
    end
end
# colvec, coef = LQOI.lqs_getrows(m, rowvec)
# TODO improve
function LQOI.lqs_getrows(lp::GLPK.Prob, idx)
    colidx, coefs = GLPK.get_mat_row(lp, idx)
    return colidx-1, coefs
end

# LQOI.lqs_getcoef(m, row, col) #??
# TODO improve
function LQOI.lqs_getcoef(m::GLPK.Prob, row, col)
    colidx, coefs = GLPK.get_mat_row(lp, row)
    idx = findfirst(colidx, col)
    if idx > 0
        return coefs[idx]
    else
        return 0.0
    end
end

# LQOI.lqs_chgcoef!(m, row, col, coef)
# TODO SPLIT THIS ONE
LQOI.lqs_chgcoef!(m::GLPK.Prob, row, col, coef)  = error("cant set singe coeff please reset full line/column")

# LQOI.lqs_delrows!(m, row, row)
function LQOI.lqs_delrows!(m::GLPK.Prob, rowbeg, rowend)

    idx = collect(rowbeg:rowend)

    GLPK.std_basis(m.inner)
    GLPK.del_rows(m.inner, length(idx), idx)

    nothing
end
# LQOI.lqs_chgctype!(m, colvec, typevec)
# TODO fix types
LQOI.lqs_chgctype!(m::GLPK.Prob, colvec, typevec) = error("changing coltype disabled")

# LQOI.lqs_chgsense!(m, rowvec, sensevec)
# TODO fix types
function LQOI.lqs_chgsense!(m::GLPK.Prob, rowvec, sensevec)
    for i in eachindex(rowvec)
        changesense!(m::GLPK.Prob, rowvec[i], sensevec[i])
    end
    nothing
end
function changesense!(m::GLPK.Prob, row, sense)
    oldsense = get_row_type(m, row)
    newsense = translatesense(sense)
    if oldsense == newsense
        return nothing
    end
    
    if newsense == GLPK.FX
        rowub = rowlb
    elseif oldsense == GLPK.DB
        if newsense == GLPK.UP
            rowlb = -Inf
        else
            rowub = Inf
        end
    else
        rowlb = get_row_lb(m, row)
        rowub = get_row_ub(m, row)
        if newsense == GLPK.UP
            rowub = rowlb
            rowlb = -Inf
        else
            rowlb = rowub
            rowub = Inf
        end
    end

    GLPK.set_row_bnds(m, row, newsense, rowlb, rowub)

    nothing 
end
function translatesense(sense)
    if sense == Cchar('E')
        return GLPK.FX
    elseif sense == Cchar('R')
        return GLPK.DB
    elseif sense == Cchar('L')
        return GLPK.UP
    elseif sense == Cchar('G')
        return GLPK.LO
    else
        error("invalid sense")
    end
end

const VAR_TYPE_MAP = Dict{Symbol,Cchar}(
    :CONTINUOUS => Cchar('C'),
    :INTEGER => Cchar('I'),
    :BINARY => Cchar('B')
)
LQOI.lqs_vartype_map(m::GLPKSolverInstance) = VAR_TYPE_MAP

# LQOI.lqs_addsos(m, colvec, valvec, typ)
LQOI.lqs_addsos!(m::GLPK.Prob, colvec, valvec, typ) = GLPK.add_sos!(m::GLPK.Prob, typ, colvec, valvec)
# LQOI.lqs_delsos(m, idx, idx)
LQOI.lqs_delsos!(m::GLPK.Prob, idx1, idx2) = GLPK.cpx_delsos!(m::GLPK.Prob, idx1, idx2)

const SOS_TYPE_MAP = Dict{Symbol,Symbol}(
    :SOS1 => :SOS1,#Cchar('1'),
    :SOS2 => :SOS2#Cchar('2')
)
LQOI.lqs_sertype_map(m::GLPKSolverInstance) = SOS_TYPE_MAP

# LQOI.lqs_getsos(m, idx)
# TODO improve getting processes
function LQOI.lqs_getsos(m::GLPK.Prob, idx)
    indices, weights, types = GLPK.cpx_getsos(m::GLPK.Prob, idx)

    return indices, weights, types == Cchar('1') ? :SOS1 : :SOS2
end
# LQOI.lqs_getnumqconstrs(m)
LQOI.lqs_getnumqconstrs(m::GLPK.Prob) = error("GLPK does not support quadratic ocnstraints")

# LQOI.lqs_addqconstr(m, cols,coefs,rhs,sense, I,J,V)
LQOI.lqs_addqconstr!(m::GLPK.Prob, cols,coefs,rhs,sense, I,J,V) = error("GLPK does not support quadratic ocnstraints")

# LQOI.lqs_chgrngval
LQOI.lqs_chgrngval!(m::GLPK.Prob, rows, vals) = GLPK.cpx_chgrngval!(m::GLPK.Prob, rows, vals)

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
LQOI.lqs_copyquad!(m::GLPK.Prob, I, J, V) = error("GLPK does no support quadratics")

# LQOI.lqs_chgobj(m, colvec,coefvec)
function LQOI.lqs_chgobj!(m::GLPK.Prob, colvec, coefvec)
    for i in eachindex(colvec)
        GLPK.set_obj_coef(m, colvec[i], coefvec[i])
    end
    nothing
end

# LQOI.lqs_chgobjsen(m, symbol)
# TODO improve min max names
function LQOI.lqs_chgobjsen!(m::GLPK.Prob, sense) 
    if sense == :Min
        GLPK.set_obj_dir(m, GLPK.MIN)
    elseif sense == :Max
        GLPK.set_obj_dir(m, GLPK.MAX)
    else
        error("Unrecognized objective sense $sense")
    end
end
# LQOI.lqs_getobj(m)
function LQOI.lqs_getobj(m::GLPK.Prob)
    n = GLPK.get_num_cols(m)
    obj = Array{Float64}(n)
    for c = 1:n
        l = GLPK.get_obj_coef(m, c)
        obj[c] = l
    end
    return obj
end

# lqs_getobjsen(m)
function LQOI.lqs_getobjsen(m::GLPK.Prob)

    s = GLPK.get_obj_dir(m)
    if s == GLPK.MIN
        return MOI.MinSense
    elseif s == GLPK.MAX
        return MOI.MaxSense
    else
        error("Internal library error")
    end
end


#=
    Variables
=#

# LQOI.lqs_getnumcols(m)
LQOI.lqs_getnumcols(m::GLPK.Prob) = GLPK.get_num_cols(m::GLPK.Prob)

# LQOI.lqs_newcols!(m, int)
LQOI.lqs_newcols!(m::GLPK.Prob, int) = GLPK.add_cols(m::GLPK.Prob, int)

# LQOI.lqs_delcols!(m, col, col)
function LQOI.lqs_delcols!(m::GLPK.Prob, col, col2)
    idx = collect(col:col2)
    GLPK.std_basis(m)
    GLPK.del_cols(m, length(idx), idx)
end

# LQOI.lqs_addmipstarts(m, colvec, valvec)
LQOI.lqs_addmipstarts!(m::GLPK.Prob, colvec, valvec) = GLPK.cpx_addmipstarts!(m::GLPK.Prob, colvec, valvec) 

#=
    Solve
=#

# LQOI.lqs_mipopt!(m)
LQOI.lqs_mipopt!(m::GLPK.Prob) = GLPK.cpx_mipopt!(m::GLPK.Prob)

# LQOI.lqs_qpopt!(m)
LQOI.lqs_qpopt!(m::GLPK.Prob) = GLPK.cpx_qpopt!(m::GLPK.Prob)

# LQOI.lqs_lpopt!(m)
LQOI.lqs_lpopt!(m::GLPK.Prob) = GLPK.cpx_lpopt!(m::GLPK.Prob)


const TERMINATION_STATUS_MAP = Dict(
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
LQOI.lqs_getx!(m::GLPK.Prob, place) = GLPK.cpx_getx!(m::GLPK.Prob, place) 

# LQOI.lqs_getax!(m, place)
LQOI.lqs_getax!(m::GLPK.Prob, place) = GLPK.cpx_getax!(m::GLPK.Prob, place)

# LQOI.lqs_getdj!(m, place)
LQOI.lqs_getdj!(m::GLPK.Prob, place) = GLPK.cpx_getdj!(m::GLPK.Prob, place)

# LQOI.lqs_getpi!(m, place)
LQOI.lqs_getpi!(m::GLPK.Prob, place) = GLPK.cpx_getpi!(m::GLPK.Prob, place)

# LQOI.lqs_getobjval(m)
LQOI.lqs_getobjval(m::GLPK.Prob) = GLPK.cpx_getobjval(m::GLPK.Prob)

# LQOI.lqs_getbestobjval(m)
LQOI.lqs_getbestobjval(m::GLPK.Prob) = GLPK.cpx_getbestobjval(m::GLPK.Prob)

# LQOI.lqs_getmiprelgap(m)
LQOI.lqs_getmiprelgap(m::GLPK.Prob) = GLPK.cpx_getmiprelgap(m::GLPK.Prob)

# LQOI.lqs_getitcnt(m)
LQOI.lqs_getitcnt(m::GLPK.Prob)  = GLPK.cpx_getitcnt(m::GLPK.Prob)

# LQOI.lqs_getbaritcnt(m)
LQOI.lqs_getbaritcnt(m::GLPK.Prob) = GLPK.cpx_getbaritcnt(m::GLPK.Prob)

# LQOI.lqs_getnodecnt(m)
LQOI.lqs_getnodecnt(m::GLPK.Prob) = GLPK.cpx_getnodecnt(m::GLPK.Prob)

# LQOI.lqs_dualfarkas(m, place)
LQOI.lqs_dualfarkas!(m::GLPK.Prob, place) = GLPK.cpx_dualfarkas!(m::GLPK.Prob, place)

# LQOI.lqs_getray(m, place)
LQOI.lqs_getray!(m::GLPK.Prob, place) = GLPK.cpx_getray!(m::GLPK.Prob, place)


MOI.free!(m::GLPKSolverInstance) = GLPK.free_model(m.inner)

"""
    writeproblem(m::AbstractSolverInstance, filename::String)
Writes the current problem data to the given file.
Supported file types are solver-dependent.
"""
writeproblem(m::GLPKSolverInstance, filename::String, flags::String="") = GLPK.write_model(m.inner, filename)


LQOI.lqs_make_problem_type_continuous(m::GLPK.Prob) = GLPK._make_problem_type_continuous(m)


#=
    old helpers
=#
function optimize!(lpm::GLPKSolverMIPInstance)
    vartype = getvartype(lpm)
    lb = getvarLB(lpm)
    ub = getvarUB(lpm)
    old_lb = copy(lb)
    old_ub = copy(ub)
    for c in 1:length(vartype)
        vartype[c] in [:Int,:Bin] && (lb[c] = ceil(lb[c]); ub[c] = floor(ub[c]))
        vartype[c] == :Bin && (lb[c] = max(lb[c],0.0); ub[c] = min(ub[c],1.0))
    end
    #lpm.cbdata.vartype = vartype
    try
        setvarLB!(lpm, lb)
        setvarUB!(lpm, ub)
        if lpm.param.presolve == GLPK.OFF
            ret_ps = GLPK.simplex(lpm.inner, lpm.smplxparam)
            ret_ps != 0 && return ret_ps
        end
        ret = GLPK.intopt(lpm.inner, lpm.param)
        if ret == GLPK.EMIPGAP || ret == GLPK.ETMLIM || ret == GLPK.ESTOP
            lpm.userlimit = true
        end
    finally
        setvarLB!(lpm, old_lb)
        setvarUB!(lpm, old_ub)
    end
end

function getvarLB(lpm::GLPKSolverInstance)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)
    lb = Array(Float64, n)
    for c = 1:n
        l = GLPK.get_col_lb(lp, c)
        if l <= -realmax(Float64)
            l = -Inf
        end
        lb[c] = l
    end
    return lb
end

function setvarLB!(lpm::GLPKSolverInstance, collb)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)
    if nonnull(collb) && length(collb) != n
        error("invalid size of collb")
    end
    for c = 1:n
        u = GLPK.get_col_ub(lp, c)
        if u >= realmax(Float64)
            u = Inf
        end
        if nonnull(collb) && collb[c] != -Inf
            l = collb[c]
            if u < Inf
                if l != u
                    GLPK.set_col_bnds(lp, c, GLPK.DB, l, u)
                else
                    GLPK.set_col_bnds(lp, c, GLPK.FX, l, u)
                end
            else
                GLPK.set_col_bnds(lp, c, GLPK.LO, l, 0.0)
            end
        else
            if u < Inf
                GLPK.set_col_bnds(lp, c, GLPK.UP, 0.0, u)
            else
                GLPK.set_col_bnds(lp, c, GLPK.FR, 0.0, 0.0)
            end
        end
    end
end

function getvarUB(lpm::GLPKSolverInstance)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)
    ub = Array(Float64, n)
    for c = 1:n
        u = GLPK.get_col_ub(lp, c)
        if u >= realmax(Float64)
            u = Inf
        end
        ub[c] = u
    end
    return ub
end

function setvarUB!(lpm::GLPKSolverInstance, colub)
    lp = lpm.inner
    n = GLPK.get_num_cols(lp)
    if nonnull(colub) && length(colub) != n
        error("invalid size of colub")
    end
    for c = 1:n
        l = GLPK.get_col_lb(lp, c)
        if l <= -realmax(Float64)
            l = -Inf
        end
        if nonnull(colub) && colub[c] != Inf
            u = colub[c]
            if l > -Inf
                if l != u
                    GLPK.set_col_bnds(lp, c, GLPK.DB, l, u)
                else
                    GLPK.set_col_bnds(lp, c, GLPK.FX, l, u)
                end
            else
                GLPK.set_col_bnds(lp, c, GLPK.UP, 0.0, u)
            end
        else
            if l > -Inf
                GLPK.set_col_bnds(lp, c, GLPK.LO, l, 0.0)
            else
                GLPK.set_col_bnds(lp, c, GLPK.FR, 0.0, 0.0)
            end
        end
    end
end






end # module