# Other possible solvers: GLPK, SCS, ECOS

function rcsMixMatrix(posRCS_B, dirRCS_B, forceRCS=ones(length(posRCS_B)), posRef_B=zeros(3))
    n = length(posRCS_B)
    dir = normalize.(dirRCS_B)
    pos = posRCS_B .- Ref(posRef_B)
    My = zeros(6, n)
    for i in 1:n
        My[:, i] = forceRCS[i]*[dir[i]; pos[i] × dir[i]]
    end
    return My
end

function rcsAllocation(u, My, c=ones(1, size(My, 2)))
    n = size(My, 2)
    y = Variable(n)
    problem = minimize(c'*y, [My*y == u, y ≤ ones(n), y ≥ zeros(n)])
    Convex.solve!(problem, ECOS.Optimizer(); silent_solver=true)
    return y.value
end

function rcsAnalysis(My; unconstrained=false)

    # Compute pure force & torque along X, Y, and Z
    n = size(My, 2)
    Umax = zeros(6); Emax = zeros(6); Umin = zeros(6); Emin = zeros(6)
    y = Variable(n)
    Mk = zeros(6, n)
    for idx in 1:6
        c = zeros(6)
        c[idx] = 1.0
        if !unconstrained
            Mk .= (I - diagm(c))*My
        end

        problem = maximize(c'*My*y, [Mk*y == zeros(6), y ≤ ones(n), y ≥ zeros(n)])
        Convex.solve!(problem, ECOS.Optimizer(); silent_solver=true)
        Umax[idx] = problem.optval
        Emax[idx] = Umax[idx]./sum(y.value.*forceRCS)

        problem = minimize(c'*My*y, [Mk*y == zeros(6), y ≤ ones(n), y ≥ zeros(n)])
        Convex.solve!(problem, ECOS.Optimizer(); silent_solver=true)
        Umin[idx] = problem.optval
        Emin[idx] = -Umin[idx]./sum(y.value.*forceRCS)
    end

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

function sampleEvenlySphere(N = 100; includeAxes=true)
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

function rcsEnvelope(My; N=1000, unconstrained=false)
    u = sampleEvenlySphere(N; includeAxes=true)

    n = size(My, 2)
    Uforce = 0.0*u
    Eforce = copy(Uforce)
    Utorque = copy(Uforce)
    Etorque = copy(Uforce)
    y = Variable(n)

    for k in eachindex(u)
        # Force envelope
        if unconstrained
            problem = maximize(u[k]'*My[1:3, :]*y, [y ≤ ones(n), y ≥ zeros(n)])
        else
            problem = maximize(u[k]'*My[1:3, :]*y, [(I - u[k]*u[k]')*My[1:3, :]*y == zeros(3), My[4:6, :]*y == zeros(3), y ≤ ones(n), y ≥ zeros(n)])
        end
        Convex.solve!(problem, ECOS.Optimizer(); silent_solver=true)
        Uforce[k] .= problem.optval*u[k]
        Eforce[k] .= Uforce[k]./sum(y.value.*forceRCS)

        # Torque envelope
        if unconstrained
            problem = maximize(u[k]'*My[4:6, :]*y, [y ≤ ones(n), y ≥ zeros(n)])
        else
            problem = maximize(u[k]'*My[4:6, :]*y, [(I - u[k]*u[k]')*My[4:6, :]*y == zeros(3), My[1:3, :]*y == zeros(3), y ≤ ones(n), y ≥ zeros(n)])
        end
        Convex.solve!(problem, ECOS.Optimizer(); silent_solver=true)
        Utorque[k] .= problem.optval*u[k]
        Etorque[k] .= Utorque[k]./sum(y.value.*forceRCS)
    end

    return Uforce, Utorque, Eforce, Etorque
end

function plotRcs(posRCS_B, dirRCS_B; scale=1.0)
    f = Figure()
    ax = LScene(f[1, 1])
    plotRcs!.(ax, posRCS_B, dirRCS_B; height=scale)
    return f, ax
end

function plotRcs!(ax, posRCS_B, dirRCS_B; scale=1.0)
    plotCone!(ax, posRCS_B, -dirRCS_B; height=scale)
end
