module MathOptInterfaceGLPK

import Base.show, Base.copy

# Standard LP interface
# importall MathProgBase.SolverInterface

using GLPK
using MathOptInterface
const MOI = MathOptInterface
using LinQuadOptInterface
const LQOI = LinQuadOptInterface

# Many functions in this modulo are adapted from GLPKMathProgInterface.jl. This is the copyright notice:
## Copyright (c) 2013: Carlo Baldassi
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
## copies of the Software, and to permit persons to whom the Software is 
## furnished to do so, subject to the following conditions:
## The above copyright notice and this permission notice shall be included in 
## all copies or substantial portions of the Software.
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.

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
    for (k,v) in s.options
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


const SUPPORTED_CONSTRAINTS_LP = [
    (LQOI.Linear, LQOI.EQ),
    (LQOI.Linear, LQOI.LE),
    (LQOI.Linear, LQOI.GE),
    (LQOI.Linear, LQOI.IV),
    (LQOI.SinVar, LQOI.EQ),
    (LQOI.SinVar, LQOI.LE),
    (LQOI.SinVar, LQOI.GE),
    (LQOI.SinVar, LQOI.IV),
    (LQOI.VecVar, MOI.Nonnegatives),
    (LQOI.VecVar, MOI.Nonpositives),
    (LQOI.VecVar, MOI.Zeros),
    (LQOI.VecLin, MOI.Nonnegatives),
    (LQOI.VecLin, MOI.Nonpositives),
    (LQOI.VecLin, MOI.Zeros)
]
const SUPPORTED_CONSTRAINTS_MIP = vcat(SUPPORTED_CONSTRAINTS_LP, [
    (LQOI.SinVar, MOI.ZeroOne),
    (LQOI.SinVar, MOI.Integer),
    # (VecVar, MOI.SOS1),
    # (VecVar, MOI.SOS2),
])
const SUPPORTED_OBJECTIVES = [
    LQOI.Linear
]

LQOI.lqs_supported_constraints(s::GLPKSolverMIP) = SUPPORTED_CONSTRAINTS_MIP
LQOI.lqs_supported_constraints(s::GLPKSolverMIPInstance) = SUPPORTED_CONSTRAINTS_MIP
LQOI.lqs_supported_objectives(s::GLPKSolverMIP) = SUPPORTED_OBJECTIVES
LQOI.lqs_supported_objectives(s::GLPKSolverMIPInstance) = SUPPORTED_OBJECTIVES

LQOI.lqs_supported_constraints(s::GLPKSolverLP) = SUPPORTED_CONSTRAINTS_LP
LQOI.lqs_supported_constraints(s::GLPKSolverLPInstance) = SUPPORTED_CONSTRAINTS_LP
LQOI.lqs_supported_objectives(s::GLPKSolverLP) = SUPPORTED_OBJECTIVES
LQOI.lqs_supported_objectives(s::GLPKSolverLPInstance) = SUPPORTED_OBJECTIVES
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
function LQOI.lqs_chgbds!(instance::GLPKSolverInstance, colvec, valvec, sensevec)
    m = instance.inner
    for i in eachindex(colvec)
        colub = Inf
        collb = -Inf
        bt = GLPK.DB
        if sensevec[i] == Cchar('E')
            colub = valvec[i]
            collb = valvec[i]
            bt = GLPK.FX
        elseif sensevec[i] == Cchar('L')
            collb = valvec[i]
            colub = Inf
            u = GLPK.get_col_ub(m, colvec[i])
            if u < Inf
                bt = GLPK.DB
                colub = u
            else
                bt = GLPK.LO
            end
        elseif sensevec[i] == Cchar('U')
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
        elseif colub == collb
            bt = GLPK.FX
        end
        GLPK.set_col_bnds(m, colvec[i], bt, collb, colub)
    end 

end


# LQOI.lqs_getlb(m, col)
LQOI.lqs_getlb(instance::GLPKSolverInstance, col) = GLPK.get_col_lb(instance.inner, col)

# LQOI.lqs_getub(m, col)
LQOI.lqs_getub(instance::GLPKSolverInstance, col) = GLPK.get_col_ub(instance.inner, col)

# LQOI.lqs_getnumrows(m)
LQOI.lqs_getnumrows(instance::GLPKSolverInstance) = GLPK.get_num_rows(instance.inner)

# LQOI.lqs_addrows!(m, rowvec, colvec, coefvec, sensevec, rhsvec)
function LQOI.lqs_addrows!(instance::GLPKSolverInstance, rowvec, colvec, coefvec, sensevec, rhsvec) 
    m = instance.inner

    nrows = length(rhsvec)

    if nrows <= 0
        error("nor row to be added")
    elseif nrows == 1
        addrow!(m, colvec::Vector, coefvec::Vector, sensevec[1], rhsvec[1])

    else

        push!(rowvec, length(colvec)+1)

        for i in 1:(length(rowvec)-1)
            inds = colvec[rowvec[i]:rowvec[i+1]-1]
            coefs = coefvec[rowvec[i]:rowvec[i+1]-1]
            addrow!(m, inds, coefs, sensevec[i], rhsvec[i])
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
    elseif sense == Cchar('R')
        # start with lower
        bt = GLPK.DB
        rowlb = rhs
        rowub = Inf
    else
        error("row type $(sense) not valid")
        bt = GLPK.FR
    end
    GLPK.set_row_bnds(lp, m, bt, rowlb, rowub)
    return
end

function setrhs(instance::GLPKSolverInstance, idx::Integer, rhs::Real)
    lp = instance.inner

    l = GLPK.get_row_lb(lp, idx)
    u = GLPK.get_row_ub(lp, idx)

    if l == u
        bt = GLPK.FX
        rowlb = rhs
        rowub = rhs
    elseif l < Inf && u < Inf
        bt = GLPK.FX
        rowlb = rhs
        rowub = rhs
    elseif l < Inf
        bt = GLPK.LO
        rowlb = rhs
        ruwub = Inf
    elseif u < Inf
        bt = GLPK.UP
        rowlb = -Inf
        ruwub = rhs
    else
        error("not valid rhs")
    end

    GLPK.set_row_bnds(lp, idx, bt, rowlb, rowub)
end
# LQOI.lqs_getrhs(m, rowvec)
function LQOI.lqs_getrhs(instance::GLPKSolverInstance, row)
    m = instance.inner
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
function LQOI.lqs_getrows(instance::GLPKSolverInstance, idx)
    lp = instance.inner
    colidx, coefs = GLPK.get_mat_row(lp, idx)
    return colidx-1, coefs
end

# LQOI.lqs_getcoef(m, row, col) #??
# TODO improve
function LQOI.lqs_getcoef(instance::GLPKSolverInstance, row, col)
    lp = instance.inner
    colidx, coefs = GLPK.get_mat_row(lp, row)
    idx = findfirst(colidx, col)
    if idx > 0
        return coefs[idx]
    else
        return 0.0
    end
end

# LQOI.lqs_chgcoef!(m, row, col, coef)
function LQOI.lqs_chgcoef!(instance::GLPKSolverInstance, row, col, coef)
    if row == 0 
        lp = instance.inner
        GLPK.set_obj_coef(lp, col, coef)
    elseif col == 0
        setrhs(instance, row, coef)
    else
        chgmatcoef!(instance, row, col, coef)
    end
    return nothing
end

function chgmatcoef!(instance::GLPKSolverInstance, row, col, coef)
    lp = instance.inner
    colidx, coefs = GLPK.get_mat_row(lp, row)
    idx = findfirst(colidx, col)
    if idx > 0
        coefs[idx] = coef
    else
        push!(colidx, col)
        push!(coefs, coef)
    end  
    GLPK.set_mat_row(lp, row, colidx, coefs)
    return nothing
end
# LQOI.lqs_delrows!(m, row, row)
function LQOI.lqs_delrows!(instance::GLPKSolverInstance, rowbeg, rowend)

    m = instance.inner

    idx = collect(rowbeg:rowend)

    GLPK.std_basis(m)
    GLPK.del_rows(m, length(idx), idx)

    nothing
end
# LQOI.lqs_chgctype!(m, colvec, typevec)
# TODO fix types
function LQOI.lqs_chgctype!(instance::GLPKSolverInstance, colvec, vartype)

    lp = instance.inner
    coltype = GLPK.CV
    for i in eachindex(colvec)
        if vartype[i] == Cint('I')
            coltype = GLPK.IV
        elseif vartype[i] == Cint('C')
            coltype = GLPK.CV
        elseif vartype[i] == Cint('B')
            coltype = GLPK.IV
            GLPK.set_col_bnds(lp, colvec[i], GLPK.DB, 0.0, 1.0)
        else
            error("invalid variable type: $(vartype[i])")
        end
        GLPK.set_col_kind(lp, colvec[i], coltype)
    end

end
# LQOI.lqs_chgsense!(m, rowvec, sensevec)
# TODO fix types
function LQOI.lqs_chgsense!(instance::GLPKSolverInstance, rowvec, sensevec)
    for i in eachindex(rowvec)
        changesense!(instance, rowvec[i], sensevec[i])
    end
    nothing
end
function changesense!(instance::GLPKSolverInstance, row, sense)
    m = instance.inner
    oldsense = GLPK.get_row_type(m, row)
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
LQOI.lqs_addsos!(instance::GLPKSolverInstance, colvec, valvec, typ) = GLPK.add_sos!(instance.inner, typ, colvec, valvec)
# LQOI.lqs_delsos(m, idx, idx)
LQOI.lqs_delsos!(instance::GLPKSolverInstance, idx1, idx2) = error("cant del SOS")

const SOS_TYPE_MAP = Dict{Symbol,Symbol}(
    :SOS1 => :SOS1,#Cchar('1'),
    :SOS2 => :SOS2#Cchar('2')
)
LQOI.lqs_sertype_map(m::GLPKSolverInstance) = SOS_TYPE_MAP

# LQOI.lqs_getsos(m, idx)
# TODO improve getting processes
function LQOI.lqs_getsos(instance::GLPKSolverInstance, idx)
    indices, weights, types = GLPK.getsos(instance.inner, idx)

    return indices, weights, types == Cchar('1') ? :SOS1 : :SOS2
end
# LQOI.lqs_getnumqconstrs(m)
LQOI.lqs_getnumqconstrs(instance::GLPKSolverInstance) = error("GLPK does not support quadratic ocnstraints")

# LQOI.lqs_addqconstr(m, cols,coefs,rhs,sense, I,J,V)
LQOI.lqs_addqconstr!(instance::GLPKSolverInstance, cols,coefs,rhs,sense, I,J,V) = error("GLPK does not support quadratic ocnstraints")

# LQOI.lqs_chgrngval
function LQOI.lqs_chgrngval!(instance::GLPKSolverInstance, rows, vals)
    lp = instance.inner
    for i in eachindex(rows)
        l = GLPK.get_row_lb(lp, rows[i])
        GLPK.set_row_bnds(lp, rows[i], GLPK.DB, l, l+vals[i])
    end
    nothing
end
const CTR_TYPE_MAP = Dict{Symbol,Cchar}(
    :RANGE => Cchar('R'),
    :LOWER => Cchar('G'),
    :UPPER => Cchar('L'),
    :EQUALITY => Cchar('E')
)
LQOI.lqs_ctrtype_map(m::GLPKSolverInstance) = CTR_TYPE_MAP

#=
    Objective
=#

# LQOI.lqs_copyquad(m, intvec,intvec, floatvec) #?
LQOI.lqs_copyquad!(instance::GLPKSolverInstance, I, J, V) = error("GLPK does no support quadratics")

# LQOI.lqs_chgobj(m, colvec,coefvec)
function LQOI.lqs_chgobj!(instance::GLPKSolverInstance, colvec, coefvec)
    m = instance.inner
    for i in eachindex(colvec)
        GLPK.set_obj_coef(m, colvec[i], coefvec[i])
    end
    nothing
end

# LQOI.lqs_chgobjsen(m, symbol)
# TODO improve min max names
function LQOI.lqs_chgobjsen!(instance::GLPKSolverInstance, sense) 
    m = instance.inner
    if sense == :Min
        GLPK.set_obj_dir(m, GLPK.MIN)
    elseif sense == :Max
        GLPK.set_obj_dir(m, GLPK.MAX)
    else
        error("Unrecognized objective sense $sense")
    end
end
# LQOI.lqs_getobj(m)
function LQOI.lqs_getobj(instance::GLPKSolverInstance)
    m = instance.inner
    n = GLPK.get_num_cols(m)
    obj = Array{Float64}(n)
    for c = 1:n
        l = GLPK.get_obj_coef(m, c)
        obj[c] = l
    end
    return obj
end

# lqs_getobjsen(m)
function LQOI.lqs_getobjsen(instance::GLPKSolverInstance)

    s = GLPK.get_obj_dir(instance.inner)
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
LQOI.lqs_getnumcols(instance::GLPKSolverInstance) = GLPK.get_num_cols(instance.inner)

# LQOI.lqs_newcols!(m, int)
function LQOI.lqs_newcols!(instance::GLPKSolverInstance, int)
    n = GLPK.get_num_cols(instance.inner)
    GLPK.add_cols(instance.inner, int)
    for i in 1:int
        GLPK.set_col_bnds(instance.inner, n+i, GLPK.FR, -Inf, Inf)
    end
    nothing
end

# LQOI.lqs_delcols!(m, col, col)
function LQOI.lqs_delcols!(instance::GLPKSolverInstance, col, col2)
    idx = collect(col:col2)
    GLPK.std_basis(instance.inner)
    GLPK.del_cols(instance.inner, length(idx), idx)
end

# LQOI.lqs_addmipstarts(m, colvec, valvec)
LQOI.lqs_addmipstarts!(instance::GLPKSolverInstance, colvec, valvec) = nothing
#=
    Solve
=#

# LQOI.lqs_mipopt!(m)
LQOI.lqs_mipopt!(instance::GLPKSolverInstance) = opt!(instance)

# LQOI.lqs_qpopt!(m)
LQOI.lqs_qpopt!(instance::GLPKSolverInstance) = error("Quadratic solving not supported")

# LQOI.lqs_lpopt!(m)
LQOI.lqs_lpopt!(instance::GLPKSolverInstance) = opt!(instance)


const TERMINATION_STATUS_MAP = Dict(
)

# LQOI.lqs_terminationstatus(m)
function LQOI.lqs_terminationstatus(model::GLPKSolverMIPInstance)

    if model.userlimit
        return MOI.OtherLimit
    end
    s = GLPK.mip_status(model.inner)
    if s == GLPK.UNDEF
        if model.param.presolve == GLPK.OFF && GLPK.get_status(model.inner) == GLPK.NOFEAS
            return MOI.InfeasibleNoResult
        else
            return MOI.OtherError
        end
    end
    if s == GLPK.OPT
        return MOI.Success
    elseif s == GLPK.INFEAS
        return MOI.InfeasibleNoResult
    elseif s == GLPK.UNBND
        return MOI.UnboundedNoResult
    elseif s == GLPK.FEAS
        return MOI.SlowProgress
    elseif s == GLPK.NOFEAS
        return MOI.OtherError
    elseif s == GLPK.UNDEF
        return MOI.OtherError
    else
        error("internal library error")
    end
end
function LQOI.lqs_terminationstatus(model::GLPKSolverLPInstance)
    s = lp_status(model)
    if s == GLPK.OPT
        return MOI.Success
    elseif s == GLPK.INFEAS
        return MOI.InfeasibleNoResult
    elseif s == GLPK.UNBND
        return OI.UnboundedNoResult
    elseif s == GLPK.FEAS
        return MOI.SlowProgress
    elseif s == GLPK.NOFEAS
        return MOI.OtherError
    elseif s == GLPK.UNDEF
        return MOI.OtherError
    else
        error("Internal library error")
    end
end

function lp_status(lpm::GLPKSolverLPInstance)
    if lpm.method == :Simplex || lpm.method == :Exact
        get_status = GLPK.get_status
    elseif lpm.method == :InteriorPoint
        get_status = GLPK.ipt_status
    else
        error("bug")
    end

    s = get_status(lpm.inner)
end

function LQOI.lqs_primalstatus(model::GLPKSolverMIPInstance)
    m = model.inner

    s = GLPK.mip_status(model.inner)

    out = MOI.UnknownResultStatus

    if s in [GLPK.OPT, GLPK.FEAS]
        out = MOI.FeasiblePoint
    end
    return out
end
function LQOI.lqs_primalstatus(model::GLPKSolverLPInstance)
    m = model.inner

    s = lp_status(model)

    out = MOI.UnknownResultStatus

    if s in [GLPK.OPT, GLPK.FEAS]
        out = MOI.FeasiblePoint
    end
    return out
end
function LQOI.lqs_dualstatus(model::GLPKSolverMIPInstance)
    return MOI.UnknownResultStatus
end
function LQOI.lqs_dualstatus(model::GLPKSolverLPInstance)
    m = model.inner
    
    s = lp_status(model)

    out = MOI.UnknownResultStatus

    if s in [GLPK.OPT]#, GLPK.FEAS]
        out = MOI.FeasiblePoint
    end
    return out
end


# LQOI.lqs_getx!(m, place)
function LQOI.lqs_getx!(instance::GLPKSolverMIPInstance, place)
    lp = instance.inner
    for c in eachindex(place)
        place[c] = GLPK.mip_col_val(lp, c)
    end
end
function LQOI.lqs_getx!(lpm::GLPKSolverLPInstance, place)
    lp = lpm.inner

    if lpm.method == :Simplex || lpm.method == :Exact
        get_col_prim = GLPK.get_col_prim
    elseif lpm.method == :InteriorPoint
        get_col_prim = GLPK.ipt_col_prim
    else
        error("bug")
    end

    for c in eachindex(place)
        place[c] = get_col_prim(lp, c)
    end
    return nothing
end

# LQOI.lqs_getax!(m, place)
function LQOI.lqs_getax!(instance::GLPKSolverMIPInstance, place)
    lp = instance.inner
    for c in eachindex(place)
        place[c] = GLPK.mip_row_val(lp, c)
    end
end
function LQOI.lqs_getax!(lpm::GLPKSolverLPInstance, place)
    lp = lpm.inner

    if lpm.method == :Simplex || lpm.method == :Exact
        get_row_prim = GLPK.get_row_prim
    elseif lpm.method == :InteriorPoint
        get_row_prim = GLPK.ipt_row_prim
    else
        error("bug")
    end

    for r in eachindex(place)
        place[r] = get_row_prim(lp, r)
    end
    return nothing
end

# LQOI.lqs_getdj!(m, place)
# no-op
function LQOI.lqs_getdj!(instance::GLPKSolverMIPInstance, place) end

function LQOI.lqs_getdj!(lpm::GLPKSolverLPInstance, place)
    lp = lpm.inner

    if lpm.method == :Simplex || lpm.method == :Exact
        get_col_dual = GLPK.get_col_dual
    elseif lpm.method == :InteriorPoint
        get_col_dual = GLPK.ipt_col_dual
    else
        error("bug")
    end

    for c in eachindex(place)
        place[c] = get_col_dual(lp, c)
    end
    return nothing
end

# LQOI.lqs_getpi!(m, place)
#no-op
function LQOI.lqs_getpi!(instance::GLPKSolverMIPInstance, place) end
function LQOI.lqs_getpi!(lpm::GLPKSolverLPInstance, place)
    lp = lpm.inner

    if lpm.method == :Simplex || lpm.method == :Exact
        get_row_dual = GLPK.get_row_dual
    elseif lpm.method == :InteriorPoint
        get_row_dual = GLPK.ipt_row_dual
    else
        error("bug")
    end

    for r in eachindex(place)
        place[r] = get_row_dual(lp, r)
    end
    return nothing
end

# LQOI.lqs_getobjval(m)
LQOI.lqs_getobjval(instance::GLPKSolverMIPInstance) = GLPK.mip_obj_val(instance.inner)

function LQOI.lqs_getobjval(lpm::GLPKSolverLPInstance)
    if lpm.method == :Simplex || lpm.method == :Exact
        get_obj_val = GLPK.get_obj_val
    elseif lpm.method == :InteriorPoint
        get_obj_val = GLPK.ipt_obj_val
    else
        error("bug")
    end
    return get_obj_val(lpm.inner)
end

# LQOI.lqs_getbestobjval(m)
LQOI.lqs_getbestobjval(instance::GLPKSolverMIPInstance) = instance.objbound

# LQOI.lqs_getmiprelgap(m)
LQOI.lqs_getmiprelgap(instance::GLPKSolverInstance) = abs(GLPK.mip_obj_val(instance.inner)-instance.objbound)/(1e-9+GLPK.mip_obj_val(instance.inner))

# LQOI.lqs_getitcnt(m)
LQOI.lqs_getitcnt(instance::GLPKSolverInstance)  = -1

# LQOI.lqs_getbaritcnt(m)
LQOI.lqs_getbaritcnt(instance::GLPKSolverInstance) = -1

# LQOI.lqs_getnodecnt(m)
LQOI.lqs_getnodecnt(instance::GLPKSolverInstance) = -1

# LQOI.lqs_dualfarkas(m, place)
LQOI.lqs_dualfarkas!(instance::GLPKSolverInstance, place) = getinfeasibilityray(instance, place)

# LQOI.lqs_getray(m, place)
LQOI.lqs_getray!(instance::GLPKSolverInstance, place) = getunboundedray(instance, place)

#no-op
function MOI.free!(instance::GLPKSolverInstance) end

"""
    writeproblem(m::AbstractSolverInstance, filename::String)
Writes the current problem data to the given file.
Supported file types are solver-dependent.
"""
writeproblem(instance::GLPKSolverInstance, filename::String, flags::String="") = GLPK.write_model(instance.inner, filename)


LQOI.lqs_make_problem_type_continuous(instance::GLPKSolverInstance) = GLPK._make_problem_type_continuous(instance.inner)


#=
    old helpers
=#
function opt!(lpm::GLPKSolverLPInstance)
    write_lp(lpm.inner, "model.lp")
    if lpm.method == :Simplex
        solve = GLPK.simplex
    elseif lpm.method == :Exact
        solve = GLPK.exact
    elseif lpm.method == :InteriorPoint
        solve = GLPK.interior
    else
        error("bug")
    end
    return solve(lpm.inner, lpm.param)
end

function opt!(lpm::GLPKSolverMIPInstance)
    write_lp(lpm.inner, "model.lp")
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
    lb = Array{Float64}(n)
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
    ub = Array{Float64}(n)
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

const vartype_map = Dict(
    GLPK.CV => :Cont,
    GLPK.IV => :Int,
    GLPK.BV => :Bin
)

function getvartype(lpm::GLPKSolverInstance)
    lp = lpm.inner
    ncol = GLPK.get_num_cols(lp)
    coltype = Array{Symbol}(ncol)
    for i in 1:ncol
        ct = GLPK.get_col_kind(lp, i)
        coltype[i] = vartype_map[ct]
        if i in lpm.binaries
            coltype[i] = :Bin
        elseif coltype[i] == :Bin # GLPK said it was binary, but we didn't tell it
            coltype[i] = :Int
        end
    end
    return coltype
end
nonnull(x) = (x != nothing && !isempty(x))

# The functions getinfeasibilityray and getunboundedray are adapted from code
# taken from the LEMON C++ optimization library. This is the copyright notice:
#
### Copyright (C) 2003-2010
### Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
### (Egervary Research Group on Combinatorial Optimization, EGRES).
###
### Permission to use, modify and distribute this software is granted
### provided that this copyright notice appears in all copies. For
### precise terms see the accompanying LICENSE file.
###
### This software is provided "AS IS" with no warranty of any kind,
### express or implied, and with no claim as to its suitability for any
### purpose.

function getinfeasibilityray(lpm::GLPKSolverLPInstance, ray)
    lp = lpm.inner

    if lpm.method == :Simplex || lpm.method == :Exact
    elseif lpm.method == :InteriorPoint
        error("getinfeasibilityray is not available when using the InteriorPoint method")
    else
        error("bug")
    end

    m = GLPK.get_num_rows(lp)

    # ray = zeros(m)
    @assert length(ray) == m

    ur = GLPK.get_unbnd_ray(lp)
    if ur != 0
        if ur <= m
            k = ur
            get_stat = GLPK.get_row_stat
            get_bind = GLPK.get_row_bind
            get_prim = GLPK.get_row_prim
            get_ub = GLPK.get_row_ub
        else
            k = ur - m
            get_stat = GLPK.get_col_stat
            get_bind = GLPK.get_col_bind
            get_prim = GLPK.get_col_prim
            get_ub = GLPK.get_col_ub
        end

        get_stat(lp, k) == GLPK.BS || error("unbounded ray is primal (use getunboundedray)")

        ray[get_bind(lp, k)] = (get_prim(lp, k) > get_ub(lp, k)) ? -1 : 1

        GLPK.btran(lp, ray)
    else
        eps = 1e-7
        for i = 1:m
            idx = GLPK.get_bhead(lp, i)
            if idx <= m
                k = idx
                get_prim = GLPK.get_row_prim
                get_ub = GLPK.get_row_ub
                get_lb = GLPK.get_row_lb
            else
                k = idx - m
                get_prim = GLPK.get_col_prim
                get_ub = GLPK.get_col_ub
                get_lb = GLPK.get_col_lb
            end

            res = get_prim(lp, k)
            if res > get_ub(lp, k) + eps
                ray[i] = -1
            elseif res < get_lb(lp, k) - eps
                ray[i] = 1
            else
                continue # ray[i] == 0
            end

            if idx <= m
                ray[i] *= GLPK.get_rii(lp, k)
            else
                ray[i] /= GLPK.get_sjj(lp, k)
            end
        end

        GLPK.btran(lp, ray)

        for i = 1:m
            ray[i] /= GLPK.get_rii(lp, i)
        end
    end

    return nothing
end

function getunboundedray(lpm::GLPKSolverLPInstance, ray)
    lp = lpm.inner

    if lpm.method == :Simplex || lpm.method == :Exact
    elseif lpm.method == :InteriorPoint
        error("getunboundedray is not available when using the InteriorPoint method")
    else
        error("bug")
    end

    m = GLPK.get_num_rows(lp)
    n = GLPK.get_num_cols(lp)

    # ray = zeros(n)
    @assert length(ray) == n

    ur = GLPK.get_unbnd_ray(lp)
    if ur != 0
        if ur <= m
            k = ur
            get_stat = GLPK.get_row_stat
            get_dual = GLPK.get_row_dual
        else
            k = ur - m
            get_stat = GLPK.get_col_stat
            get_dual = GLPK.get_col_dual
            ray[k] = 1
        end

        get_stat(lp, k) != GLPK.BS || error("unbounded ray is dual (use getinfeasibilityray)")

        for (ri, rv) in zip(GLPK.eval_tab_col(lp, ur)...)
            ri > m && (ray[ri - m] = rv)
        end

        if (GLPK.get_obj_dir(lp) == GLPK.MAX) $ (get_dual(lp, k) > 0)
            scale!(ray, -1.0)
        end
    else
        for i = 1:n
            ray[i] = GLPK.get_col_prim(lp, i)
        end
    end

    return nothing
end





end # module