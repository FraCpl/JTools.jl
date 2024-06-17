# Other possible solvers: GLPK, SCS, ECOS

"""
    My = rcsMixMatrix(posRCS, dirRCS, forceRCS, posRef)

Compute the RCS mixing matrix ```My``` with respect to ```posRef```, so that:

[Force; Torque] = ```My```*y

with 0 ≤ y ≤ 1.
"""
function rcsMixMatrix(posRCS, dirRCS, forceRCS=ones(length(posRCS)), posRef=zeros(3))
    n = length(posRCS)
    dir = normalize.(dirRCS)
    pos = posRCS .- Ref(posRef)
    My = zeros(6, n)
    for i in 1:n
        My[:, i] = forceRCS[i]*[dir[i]; pos[i] × dir[i]]
    end
    return My
end

function rcsAnalysis(My; unconstrained=false)

    v = [[+1.0; 0.0; 0.0], [0.0; +1.0; 0.0], [0.0; 0.0; +1.0],
         [-1.0; 0.0; 0.0], [0.0; -1.0; 0.0], [0.0; 0.0; -1.0]]

    Uforce, Utorque, Eforce, Etorque = rcsEnvelope(My; unconstrained=unconstrained, v=v)
    Umax = [norm.(Uforce[1:3]); norm.(Utorque[1:3])]
    Umin = -[norm.(Uforce[4:6]); norm.(Utorque[4:6])]
    Emax = [norm.(Eforce[1:3]); norm.(Etorque[1:3])]
    Emin = [norm.(Eforce[4:6]); norm.(Etorque[4:6])]

    type = unconstrained ? "Unconstrained " : "Pure "
    println(type*"Force")
    println("  X:  [$(round(Umin[1]; digits=3)), $(round(Umax[1]; digits=3))] N")
    println("  Y:  [$(round(Umin[2]; digits=3)), $(round(Umax[2]; digits=3))] N")
    println("  Z:  [$(round(Umin[3]; digits=3)), $(round(Umax[3]; digits=3))] N")
    println("  ηX: [$(round(Emin[1]; digits=3)), $(round(Emax[1]; digits=3))]")
    println("  ηY: [$(round(Emin[2]; digits=3)), $(round(Emax[2]; digits=3))]")
    println("  ηZ: [$(round(Emin[3]; digits=3)), $(round(Emax[3]; digits=3))]\n")

    println(type*"Torque")
    println("  X:  [$(round(Umin[4]; digits=3)), $(round(Umax[4]; digits=3))] Nm")
    println("  Y:  [$(round(Umin[5]; digits=3)), $(round(Umax[5]; digits=3))] Nm")
    println("  Z:  [$(round(Umin[6]; digits=3)), $(round(Umax[6]; digits=3))] Nm")
    println("  rX: [$(round(Emin[4]; digits=3)), $(round(Emax[4]; digits=3))] m")
    println("  rY: [$(round(Emin[5]; digits=3)), $(round(Emax[5]; digits=3))] m")
    println("  rZ: [$(round(Emin[6]; digits=3)), $(round(Emax[6]; digits=3))] m")
    return nothing
end

function sampleEvenlySphere(N=100; includeAxes=true)
    a = 4π/N
    d = √a
    Mθ = round(Int, π/d)
    dθ = π/Mθ
    dϕ = a/dθ

    r = []
    for m in 0:Mθ - 1
        θ = π*(m + 0.5)/Mθ
        Mϕ = round(Int, 2π*sin(θ)/dϕ)
        for n in 0:Mϕ - 1
            ϕ = 2π*n/Mϕ
            sθ, cθ = sincos(θ)
            sϕ, cϕ = sincos(ϕ)
            push!(r, [sθ*cϕ; sθ*sϕ; cθ])
        end
    end
    if includeAxes
        for i in 1:3
            u = zeros(3)
            u[i] = 1.0
            push!(r, u)
            push!(r, -u)
        end
    end
    return r
end

function rcsEnvelope(My; N=1000, unconstrained=false, v=sampleEvenlySphere(N; includeAxes=true))

    v .= normalize.(v)

    forceRCS = [norm(f[1:3]) for f in eachcol(My)]
    n = size(My, 2)
    Uforce = deepcopy(v)
    Eforce = deepcopy(v)
    Utorque = deepcopy(v)
    Etorque = deepcopy(v)
    y = Variable(n)

    for k in eachindex(v)
        # Force envelope - Maximize force along given direction v[k]
        if unconstrained
            problem = maximize(v[k]'*My[1:3, :]*y, [y ≤ ones(n), y ≥ zeros(n)])
        else
            problem = maximize(v[k]'*My[1:3, :]*y, [(I - v[k]*v[k]')*My[1:3, :]*y == zeros(3), My[4:6, :]*y == zeros(3), y ≤ ones(n), y ≥ zeros(n)])
        end
        Convex.solve!(problem, () -> ECOS.Optimizer(); silent=true)
        Uforce[k] .*= problem.optval
        Eforce[k] .*= problem.optval/sum(y.value.*forceRCS)

        # Torque envelope - Maximize torque along given direction v[k]
        if unconstrained
            problem = maximize(v[k]'*My[4:6, :]*y, [y ≤ ones(n), y ≥ zeros(n)])
        else
            problem = maximize(v[k]'*My[4:6, :]*y, [(I - v[k]*v[k]')*My[4:6, :]*y == zeros(3), My[1:3, :]*y == zeros(3), y ≤ ones(n), y ≥ zeros(n)])
        end
        Convex.solve!(problem, () -> ECOS.Optimizer(); silent=true)
        Utorque[k] .*= problem.optval
        Etorque[k] .*= problem.optval/sum(y.value.*forceRCS)
    end

    return Uforce, Utorque, Eforce, Etorque
end

function plotRcs(posRCS, dirRCS; scale=1.0)
    f = Figure(); display(f)
    ax = LScene(f[1, 1])
    plotRcs!(ax, posRCS, dirRCS; scale=scale)
    return f, ax
end

function plotRcs!(ax, posRCS, dirRCS; scale=1.0)
    plotCone!.(ax, posRCS, -dirRCS; height=scale)
end


function rcsAllocation(u, My, c=ones(size(My, 2)))
    n = size(My, 2)
    y = Variable(n)
    problem = minimize(c'*y, [My*y == u, y ≤ ones(n), y ≥ zeros(n)])
    Convex.solve!(problem, () -> ECOS.Optimizer(); silent=true)
    return y.value
end

#runTMFSimplexMIT
#	This function solves the thrusters selection problem using the
#	Simplex algorithm.
#
#   min  c'*y
#   s.t. My*y = u
#        ylb <= y <= yub
#
#   Inputs:
#   u       Torsor to be realized [Force; Torque]
#   My      Dimensional mixing matrix, i.e. u = My*y, where y in [ylb, yub]
#   c       Cost vector, typically c = F./isp (default: 1)
#   ylb     [OPT] Lower bound on y (default: 0)
#   yub     [OPT] Upper bound on y (default: 1)
#
#   Outputs:
#   y       Solution to the LP problem (i.e., thrusters opening ratios)
#   iter    Number of iterations performed
#   flag    1: solution found, 0: solution not found
#
#	References:
#	[1] T. Yang, "Optimal Thruster Selection with Robust Estimation for
#       Formation Flying Applications". PhD Thesis, MIT
#   [2] J. Paradiso, "A Highly Adaptable Steering/Selection Procedure for
#       Combined CMG/RCS Spacecraft Control".
#
#	Authors: F. Capolupo
function rcsAllocationSimplex(u, My, c=ones(size(My, 2)))
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Justification of cost 'z' and matrix 'E'.
    # If we split  the 'y' vector into basis and non-basis components, we can
    # rewrite the LP problem as:
    # J = cb*yb + cn*yn     (1.a)
    # Myb*yb + Myn*yn = u   (1.b)
    #
    # Solving for yb in (1.b) leads to
    # yb = inv(Myb)*u - inv(Myb)*Myn*yn   (2)
    #
    # Substituting (2) into (1.a)
    # J = cb*inv(Myb)*u - cb*inv(Myb)*Myn*yn + cn*yn
    #
    # The first term is useles because it does not depend from y, therefore the
    # new cost function to be minimized is
    # J' = (cn - cb*inv(Myb)*Myn)*yn
    #    = (cn - cb*E)*yn
    #
    # so that we can define z = cb*E - cn, to be maximized.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    xsign(u) = u < 0.0 ? -1.0 : 1.0

    m, n = size(My)           # m = number of dof, n = number of thrusters

    # if nargin < 6
    #     maxIter = 3*m+1;
    #     if nargin < 5
    #         yub = ones(n,1);
    #         if nargin < 4
    #             ylb = zeros(n,1);
    #             if nargin < 3
    #                 c = ones(1,n);
    #             end
    #         end
    #     end
    # end
    #
    # # Variable change to have ylb = 0
    # off = false;
    # if any(ylb ~= 0)
    #     u = u - My*ylb;
    #     yub = yub - ylb;
    #     off = true;
    # end

    maxIter = 3m + 1
    K = 1e3*maximum(c)*ones(1, m)       # Slack variables cost

    # Setup the initial solution
    Ymax = [ones(n); 1e5*ones(m)]       # 0 <= Y <= Ymax, where Y = [y; s]
    E = diagm(xsign.(u))*My 	        # [m x n] This is E = inv(My[basis])*My[nonBasis], [1]: E <-- Y = inv(B)*A
    z  = (K*E - c')[:]                  # [n] Cost change when bringing in the base a thruster which is out of the basis (i.e., increasing Yn[i] from 0)
    ib = Integer.(n+1:n+m)              # [m] Index of thrusters in the basis. Initial solution is y = s
    Yn = zeros(n + m)                   # [n+m] Thrusters out of the basis, either at zero (Yn[i] = 0) or at max (Yn[i] = Ymax[i])
    Yb = abs.(u)                        # [m] Basis vector (i.e., y of the m thrusters that form the basis)

    # Loop until the evauators are non-positive
    case2 = false
    for _ in 1:maxIter

        # Find the maximum evaluator
        # The thrusters wich maximize 'z' is invited in the basis
        # (and if zMax > 0, i.e. if there is still room for improvement)
        zMax, iMax = findmax(first, z)
        if zMax ≤ 0.0; break; end

        # Determine the candidate thrusters to leave the basis
        # (best pivot and best upper bnd index)
        e = E[:, iMax]      # Local vector indicating the rate of change of the basic thrusters due to the value of the invited thruster
        rPivot = Inf        # Smallest pivot ratio
        rUpper = Inf        # Smallest upper bound ratio
        jPivot = 0          # Local (i.e., within basis) index of the smallest pivot ratio
        jUpper = 0          # Local (i.e., within basis) index of the smallest upper bound ratio

        for j = 1:m
            if e[j] > 0
                # The thruster in the basis will go to zero and leave the basis
                # (CASE 1)
                r = Yb[j]/e[j]  # [CHECK] removed 'abs': abs(Yb(i))/e(i);
                if r < rPivot
                    rPivot = r
                    jPivot = j
                end
            elseif e[j] < 0
                # The thruster in the basis will go to its upper bound and
                # leave the basis (CASE 2)
                r = (Yb[j] - Ymax[ib[j]])/e[j]
                if r < rUpper
                    rUpper = r
                    jUpper = j
                end
            end
        end

        # Now determine how the incoming variable should be incorporated and
        # how the basic variables should be removed from basis
        if rUpper ≥ Ymax[iMax] && rPivot ≥ Ymax[iMax]
            # CASE 3: the non-basic variable stays out of the basis but goes to
            # its upper bnd. The basis is not modified.
            Yn[iMax] = Ymax[iMax]   # -Yn(iMax);      # update the non-basis vector list. Error was here: the case 1 -> 0 was not treated [CHECK: never needed?] (old code: Yout(iMax) = Ymax(iMax))
            Yb .-= Ymax[iMax]*e                 # update the on times for the jets in the basis (substract the effects of having an out of basis thruster at its upper bound)

            # Switching polarity to indicate that now thruster ratio y must be
            # decreased rather than increased when trying to decrease cost in
            # next iterations (increase z)
            E[:, iMax] .*= -1.0
            z[iMax] *= -1.0

        else
            # Update list of thrusters out of the basis
            if rUpper ≥ rPivot && rPivot ≤ Ymax[iMax]
                # CASE 1: one basic variable goes to zero and leaves the basis
                # while a non-basis variable gets added
                jOut = jPivot           # Local index of the thruster removed from the basis
                iOut = ib[jOut]   	    # Global index of the thruster removed from the basis
                r = rPivot
                Yn[iOut] = 0.0          # Thruster goes out of the basis and is shut down

            else #elseif (rPivot > rUpper) && (rUpper <= Ymax(iMax)) # [CHECK] was elseif
                # CASE2: the non-basic variable enters the basis while the one
                # it replaces goes to its upper bnd
                jOut = jUpper           # Local index of the thruster removed from the basis
                iOut = ib[jOut]         # Global index of the thruster removed from the basis
                r = rUpper
                Yn[iOut] = Ymax[iOut]   # Thruster goes out of the basis and to its upper bound
                case2 = true
            end

            # Update basis
            Yb = Yb - r*e               # Update the solution (substract the effects of having removed one thruster from  basis)
            ib[jOut] = iMax             # Update list of thrusters in basis (replace out with max)

            # Chech to see if incoming var is decreasing from upper bnd
            if Yn[iMax] != 0.0
                Yb[jOut] = Yn[iMax] - r
                Yn[iMax] = 0.0

                # Switching polarity back to 'nominal' to indicate that now
                # thruster ratio y must be increased rather than decreased from
                # its upper bound
                E[:, iMax] .*= -1.0
                e .*= -1.0
                z[iMax] *= -1.0
                zMax *= -1.0
            else
                Yb[jOut] = r
            end

            # Transform the linear combination coefficients and evaluators.
            # This corresponds to Eq. (11) of [1].
            eNew = E[jOut, :]/e[jOut]
            z .-= zMax*eNew
            E .-= e*eNew'
            E[jOut, :] = eNew

            if case2
                # Switching polarity to indicate that now thruster ratio y must
                # be decreased from its upper bound rather than increased from
                # zero.
                if iOut ≤ n
                    E[:, iOut] .*= -1.0
                    z[iOut] *= -1.0
                end
                case2 = false  # Reset case2 flag
            end
        end
    end

    # Build the final solution Y. Y includes the 'm' basic variables, plus any
    # other non-basic variable at its upper bound.
    Yn[ib] .= Yb   # Merge non-basis and basis data
    y = Yn[1:n]   # Remove slack variables from solution

    # # Restore ylb offset
    # if off, y = y + ylb; end

    # Generate success flag - All slack variables must be zero
    #flag = norm(Yn[n+1:end]) < 1e-7;
    return y
end
