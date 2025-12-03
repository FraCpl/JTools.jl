using Delaunay
using GeometryBasics
using LinearAlgebra
using JTools
using GLMakie

function analyzeRcsConf(posRCS, dirRCS, forceRCS=ones(length(posRCS)), posRef=zeros(3); N=600, unconstrained=false)
    function points2mesh(x)
        m = delaunay(mapreduce(permutedims, vcat, x))
        f = [TriangleFace([r[1], r[2], r[3]]) for r in eachrow(m.convex_hull)]

        # Check normals to fix shading
        for i in eachindex(f)
            i1, i2, i3 = f[i]
            v1 = m.points[i1, :]
            v2 = m.points[i2, :]
            v3 = m.points[i3, :]
            N = cross(v2 - v1, v3 - v2)
            d = v1 + v2 + v3
            if dot(N, d) < 0
                f[i] = [i3, i2, i1]
            end
        end
        return GeometryBasics.Mesh(Makie.to_vertices(m.points), f)
    end

    # Analysis
    My = rcsMixMatrix(posRCS, dirRCS, forceRCS, posRef)
    Umax, Umin, Emax, Emin = rcsAnalysis(My; unconstrained=unconstrained)
    Uforce, Utorque, Eforce, Etorque = rcsEnvelope(My; N=N, unconstrained=unconstrained)

    # Plot
    fig = Figure(; size=(1000, 600));
    display(fig)
    axf = Axis3(fig[1, 1]; aspect=:data, title="Pure force and efficiency (T = 1 N)")
    axt = Axis3(fig[1, 2]; aspect=:data, title="Pure torque and efficiency (T = 1 N)")
    # scatter3!(axf, Uforce; markersize=5)
    # scatter3!(axf, Eforce; markersize=5)
    # scatter3!(axt, Utorque; markersize=5)
    # scatter3!(axt, Etorque; markersize=5)

    GLMakie.mesh!(axf, points2mesh(Uforce); transparency=true, overdraw=false, color=:cyan, alpha=0.6)
    GLMakie.mesh!(axf, points2mesh(Eforce); transparency=true, overdraw=false, color=:rosybrown2, alpha=0.6)
    GLMakie.mesh!(axt, points2mesh(Utorque); transparency=true, overdraw=false, color=:cyan, alpha=0.6)
    GLMakie.mesh!(axt, points2mesh(Etorque); transparency=true, overdraw=false, color=:rosybrown2, alpha=0.6)

    titfx = "X: [$(round(Umin[1]; digits=3)), $(round(Umax[1]; digits=3))] N,    ηX: [$(round(Emin[1]; digits=3)), $(round(Emax[1]; digits=3))]"
    titfy = "Y: [$(round(Umin[2]; digits=3)), $(round(Umax[2]; digits=3))] N,    ηY: [$(round(Emin[2]; digits=3)), $(round(Emax[2]; digits=3))]"
    titfz = "Z: [$(round(Umin[3]; digits=3)), $(round(Umax[3]; digits=3))] N,    ηZ: [$(round(Emin[3]; digits=3)), $(round(Emax[3]; digits=3))]"
    tittx = "X: [$(round(Umin[4]; digits=3)), $(round(Umax[4]; digits=3))] Nm,    rX: [$(round(Emin[4]; digits=3)), $(round(Emax[4]; digits=3))]"
    titty = "Y: [$(round(Umin[5]; digits=3)), $(round(Umax[5]; digits=3))] Nm,    rY: [$(round(Emin[5]; digits=3)), $(round(Emax[5]; digits=3))]"
    tittz = "Z: [$(round(Umin[6]; digits=3)), $(round(Umax[6]; digits=3))] Nm,    rZ: [$(round(Emin[6]; digits=3)), $(round(Emax[6]; digits=3))]"
    Label(fig[2, 1], "$titfx\n$titfy\n$titfz"; tellwidth=false, halign=:left, justification=:left)
    Label(fig[2, 2], "$tittx\n$titty\n$tittz"; tellwidth=false, halign=:left, justification=:left)

    return fig
end
