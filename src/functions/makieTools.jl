function theme_fra(spinevis=true; light=false)
    if light
        bg = :white
        fg = RGBf(2/255,50/255,71/255)
        gridcol = (:black, 0.25)
        minorgridcol = (:black, 0.2)
        fga = (RGBf(2/255,50/255,71/255), 0.5)
    else
        bg = RGBf(40/255, 44/255, 52/255)
        fg = RGBf(0.7,0.7,0.7)
        gridcol = (:white, 0.09)
        minorgridcol = (:white, 0.02)
        fga = RGBf(0.3,0.3,0.3)
    end
    Theme(
        backgroundcolor = bg,
        textcolor = fg,
        linecolor = fga,
        Axis = (
            backgroundcolor = :transparent,
            xgridcolor = gridcol,
            ygridcolor = gridcol,
            leftspinevisible = spinevis,
            rightspinevisible = spinevis,
            bottomspinevisible = spinevis,
            topspinevisible = spinevis,
            leftspinecolor = fga,
            rightspinecolor = fga,
            bottomspinecolor = fga,
            topspinecolor = fga,
            xminorticksvisible = false,
            yminorticksvisible = false,
            xminorgridcolor = minorgridcol,
            yminorgridcolor = minorgridcol,
            xticksvisible = false,
            yticksvisible = false,
            xlabelpadding = 3,
            ylabelpadding = 3
        ),
        Legend = (
            framevisible = true,
            backgroundcolor=bg,
            framewidth=3.0,
            framecolor=bg,
            padding = (0, 0, 0, 0),
        ),
        Axis3 = (
            xgridcolor = gridcol,
            ygridcolor = gridcol,
            zgridcolor = gridcol,
            xspinesvisible = spinevis,
            yspinesvisible = spinevis,
            zspinesvisible = spinevis,
            xspinecolor_1 = fga,
            xspinecolor_2 = fga,
            xspinecolor_3 = fga,
            yspinecolor_1 = fga,
            yspinecolor_2 = fga,
            yspinecolor_3 = fga,
            zspinecolor_1 = fga,
            zspinecolor_2 = fga,
            zspinecolor_3 = fga,
            xticksvisible = false,
            yticksvisible = false,
            zticksvisible = false,
        ),
        Colorbar = (
            ticksvisible = false,
            spinewidth = 0,
            ticklabelpad = 5,
        )
    )
end

# Example:
#using GLMakie
#using JTools
#f = Figure()
#ax = LScene(f[1, 1]; show_axis=false)
#plotMoon!(ax)
function plotMoon!(ax, radius=1.0, pos_I=zeros(3), R_IM=I; n=256, kwargs...)
    plotPlanet!(ax, joinpath(@__DIR__, "../../assets/lroc_color_poles_2k.tif"), radius, pos_I, R_IM; n=n, kwargs...)
end

function plotEarth!(ax, radius=1.0, pos_I=zeros(3), R_IE=I; n=256, kwargs...)
    plotPlanet!(ax, joinpath(@__DIR__, "../../assets/earth.jpg"), radius, pos_I, R_IE; n=n, kwargs...)
end

function plotPlanet!(ax, texture, radius=1.0, pos_I=zeros(3), R_IP=I; n=256, kwargs...)
    θ = LinRange(0, π, n)       # Colatitude
    φ = LinRange(-π, π, 2n)     # Longitude
    posBody_I = [R_IP*(radius*[cos(φ)*sin(θ); sin(φ)*sin(θ); cos(θ)]) + pos_I  for θ in θ, φ in φ]
    surface!(ax, getindex.(posBody_I, 1), getindex.(posBody_I, 2), getindex.(posBody_I, 3);
        color=load(texture), shading=NoShading, kwargs...)
end

function plotSphere!(ax, radius=1.0, pos_I=zeros(3), n=256; kwargs...)
    θ = LinRange(0, π, n)       # Colatitude
    φ = LinRange(-π, π, 2n)     # Longitude
    posBody_I = [radius*[cos(φ)*sin(θ); sin(φ)*sin(θ); cos(θ)] + pos_I  for θ in θ, φ in φ]
    surface!(ax, getindex.(posBody_I, 1), getindex.(posBody_I, 2), getindex.(posBody_I, 3);
        shading=NoShading, kwargs...)
end

function axisoff!(ax)
    hidedecorations!(ax)  # hides ticks, grid and lables
    hidespines!(ax)  # hide the frame
end

# Plot axes of frame F in image (world) frame I
# CAUTION: This seems to make MAKIE crash
function plotframe!(ax, pos_I=zeros(3), R_IF=Matrix(1.0I, 3, 3), length=1;
        colors=[:red, :green, :blue], sub="",
        labels=[rich("x", subscript(sub)), rich("y", subscript(sub)), rich("z", subscript(sub))],
        labelspace=1.1, kwargs...)
    for i in 1:3
        u = R_IF[:, i]*length
        lines3!(ax, [pos_I, pos_I + u]; color=colors[i], kwargs...)
        text!(ax, pos_I[1] + labelspace*u[1], pos_I[2] + labelspace*u[2], pos_I[3] + labelspace*u[3]; text=labels[i], align=(:center, :center), color=colors[i])
    end
end

# Transform 3D model
function transformModel(m, pos_I, R_IB=I, scale=1.0)
    c, f = coordinates(m), faces(m)
    mt = [Point3f(pos_I + scale*R_IB*c[k]) for k in eachindex(c)]
    return GeometryBasics.Mesh(mt, f)
end

# FOV is full field of view
function sensorFovModel(; pos_I=zeros(3), R_IS=I, FOV=(40π/180, 40π/180), length=1.0)
    x = length*tan(FOV[1]/2)
    y = length*tan(FOV[2]/2)
    v = [Point3f(pos_I + [0.0, 0.0, 0.0]),
          Point3f(pos_I + R_IS*[-x, +y, length]),
          Point3f(pos_I + R_IS*[-x, -y, length]),
          Point3f(pos_I + R_IS*[+x, -y, length]),
          Point3f(pos_I + R_IS*[+x, +y, length]),]
    f = [TriangleFace([1, 2, 3]), TriangleFace([1, 3, 4]), TriangleFace([1, 4, 5]), TriangleFace([1, 5, 2])]
    return GeometryBasics.Mesh(v, f)
end

# FOV is full field of view
function plotSensorFov!(ax, pos_I=zeros(3), R_IS=I, FOV=(40π/180, 40π/180); length=1.0, kwargs...)
    mesh!(ax, sensorFovModel(pos_I=pos_I, R_IS=R_IS, FOV=FOV, length=length); kwargs...)
end

function plotCuboid!(ax, lx=1.0, ly=1.0, lz=1.0, pos_I=zeros(3), R_IB=I; kwargs...)
    m = cuboidModel(lx, ly, lz; pos_I=pos_I, R_IB=R_IB)
    mesh!(ax, m; kwargs...)
end

function cuboidModel(lx, ly, lz; pos_I=zeros(3), R_IB=I, usequad=false)
    v = [Point3f(pos_I + R_IB*[+lx/2; -ly/2; -lz/2]),
         Point3f(pos_I + R_IB*[+lx/2; +ly/2; -lz/2]),
         Point3f(pos_I + R_IB*[+lx/2; +ly/2; +lz/2]),
         Point3f(pos_I + R_IB*[+lx/2; -ly/2; +lz/2]),
         Point3f(pos_I + R_IB*[-lx/2; -ly/2; -lz/2]),
         Point3f(pos_I + R_IB*[-lx/2; +ly/2; -lz/2]),
         Point3f(pos_I + R_IB*[-lx/2; +ly/2; +lz/2]),
         Point3f(pos_I + R_IB*[-lx/2; -ly/2; +lz/2]),
    ]
    if usequad
        # use quad faces instead of triangles
        f = [
            QuadFace([1, 2, 3, 4]),
            QuadFace([2, 3, 7, 6]),
            QuadFace([7, 8, 5, 6]),
            QuadFace([1, 4, 8, 5]),
            QuadFace([4, 3, 7, 8]),
            QuadFace([1, 2, 6, 5]),
        ]
    else
        f = [
            TriangleFace([1, 2, 3]), TriangleFace([1, 3, 4]),
            TriangleFace([2, 3, 7]), TriangleFace([2, 7, 6]),
            TriangleFace([7, 8, 5]), TriangleFace([7, 5, 6]),
            TriangleFace([1, 4, 8]), TriangleFace([1, 8, 5]),
            TriangleFace([4, 3, 7]), TriangleFace([4, 7, 8]),
            TriangleFace([1, 2, 6]), TriangleFace([1, 6, 5]),
        ]
    end
    return GeometryBasics.Mesh(v, f)
end

function plotCone!(ax, posTip=zeros(3), dirAxis=[0.0; 0.0; 1.0]; height=1.0, fullAngle=40π/180, kwargs...)
    N = 200
    z = normalize(dirAxis)
    x = normalize(randn(3) × z)
    y = z × x
    R = [x y z]
    r = height*tan(fullAngle/2)
    circ = Ref(posTip) .+ Ref(R).*[[r*sin(θ); r*cos(θ); height] for θ in range(0, 2π, N)]
    lower = fill(Point3f(posTip[1], posTip[2], posTip[3]), N)
    upper = [Point3f(c[1], c[2], c[3]) for c in circ]
    band!(ax, lower, upper; kwargs...)
    lines!(ax, upper; color=:black)
end

function plotCylinder!(ax, posCenter=zeros(3), dirAxis=[0.0; 0.0; 1.0]; height=1.0, radius=1.0, noedges=false, kwargs...)
    N = 200
    z = normalize(dirAxis)
    x = normalize(randn(3) × z)
    y = z × x
    R = [x y z]

    circup = Ref(posCenter) .+ Ref(R).*[[radius*sin(θ); radius*cos(θ); +height/2] for θ in range(0, 2π, N)]
    circdw = Ref(posCenter) .+ Ref(R).*[[radius*sin(θ); radius*cos(θ); -height/2] for θ in range(0, 2π, N)]
    upper = [Point3f(c[1], c[2], c[3]) for c in circup]
    lower = [Point3f(c[1], c[2], c[3]) for c in circdw]
    band!(ax, lower, upper; kwargs...)

    circupc = posCenter + R*[0.0; 0.0; height/2]
    upc = fill(Point3f(circupc[1], circupc[2], circupc[3]), N)
    band!(ax, upc, upper; kwargs...)

    circdwn = posCenter - R*[0.0; 0.0; height/2]
    upd = fill(Point3f(circdwn[1], circdwn[2], circdwn[3]), N)
    band!(ax, upd, lower; kwargs...)

    if !noedges
        lines!(ax, upper; color=:black)
        lines!(ax, lower; color=:black)
    end
end

function plotEcdf(x, F; kwargs...)
    f = Figure()
    ax = Axis(f[1, 1]; ylabel="CDF", limits=(nothing, nothing, 0.0, nothing))
    plotEcdf!(ax, x, F; kwargs...)
    return f, ax
end

function plotEcdf!(ax, x, F; kwargs...)
    dd = hist!(ax, x; bins=minimum([length(unique(x)); 150]), normalization=:pdf, scale_to=1.0, alpha=0.3, strokewidth=0.0, kwargs...)
    dd.label = nothing
    stairs!(x, F; kwargs...)
    return nothing
end

function plotBox(i, x, width=0.3)
    f = Figure()
    ax = Axis(f[1, 1])
    plotBox!(ax, i, x; width=width)
    return f, ax
end

function plotBox!(ax, i, x; width=0.3, kwargs...)
    xF, F = ecdf(x)

    μ = xF[findfirst(F .≥ 0.5)]# sum(x)/length(x)
    x1up = xF[findfirst(F .≥ 0.75)]
    x1dw = xF[findfirst(F .≥ 0.25)]
    IR = x1up - x1dw
    x3up = x1up + 1.5IR # xF[findfirst(F .≥ 0.997)]
    x3dw = x1dw - 1.5IR # xF[findfirst(F .≥ 1 - 0.997)]

    violin!(ax, 0*x .+ i, x; datalimits=extrema, kwargs...)
    #boxplot!(ax,0*x .+ i, x; width=0.5)

    lines!(ax, [-1.0; 1.0]*width/2 .+ i, [μ; μ]; color=:black, linewidth=2)
    lines!(ax, [-1.0; 1.0; 1.0; -1.0; -1.0]*width/2 .+ i, [x1dw; x1dw; x1up; x1up; x1dw]; color=:black, linewidth=1)
    lines!(ax, [i; i], [x1up; x3up]; color=:black, linewidth=1.5)
    lines!(ax, [i; i], [x1dw; x3dw]; color=:black, linewidth=1.5)
    lines!(ax, [-1.0; 1.0]*width/7 .+ i, [x3up; x3up]; color=:black, linewidth=2)
    lines!(ax, [-1.0; 1.0]*width/7 .+ i, [x3dw; x3dw]; color=:black, linewidth=2)

    #p_big = decompose(Point2f, Circle(Point2f(0), 1))
    #p_small = decompose(Point2f, Circle(Point2f(0), 0.5))
    #marker = Polygon(p_big, [p_small])
    iout = findall(x .> x3up .|| x.< x3dw)
    scatter!(ax, 0*x[iout] .+ i, x[iout]; marker=:cross, color=:black, markersize=5)

    return nothing
end

function lines3!(ax, pos; kwargs...)
    lines!(ax, getindex.(pos, 1), getindex.(pos, 2), getindex.(pos, 3); kwargs...)
end

function scatter3!(ax, pos; kwargs...)
    scatter!(ax, getindex.(pos, 1), getindex.(pos, 2), getindex.(pos, 3); kwargs...)
end

function multilines!(ax, x, y, idx=1:lastindex(y[1]); colors=0.0, kwargs...)
    for i in idx
        if colors == 0.0
            lines!(ax, x, getindex.(y, i); kwargs...)
        else
            lines!(ax, x, getindex.(y, i); color=colors[i], kwargs...)
        end
    end
end

function mergeMesh(m1, m2)
    f = faces(m1)
    v = [m1.position; m2.position]
    nv = length(m1.position)
    for f2 in faces(m2)
        push!(f, TriangleFace([f2[1]+nv, f2[2]+nv, f2[3]+nv]))
    end
    return GeometryBasics.Mesh(v, f)
end

# Example input data:
# data = (A = randn(n), B = randn(n), C = randn(n), D = randn(n))
# TODO: Synchronize zoom among subplots
function plotCorrelation(data, kwargs...)
    fig = Figure(size=(800, 800)); display(fig)
    plotCorrelation!(fig, data, kwargs...)
    return fig
end

function plotCorrelation!(fig, data, kwargs...)

    varnames = keys(data)
    nvars = length(varnames)

    # Loop over variable pairs to plot
    for i in 1:nvars, j in 1:nvars
        if i > j
            if j == 1
                if i == nvars
                    ax = Axis(fig[i-1, j]; ylabel="$(varnames[i])", xlabel="$(varnames[j])")
                else
                    ax = Axis(fig[i-1, j]; ylabel="$(varnames[i])", xticksvisible=false, xticklabelsvisible=false)#
                end
            elseif i == nvars
                ax = Axis(fig[i-1, j]; xlabel="$(varnames[j])", yticksvisible=false, yticklabelsvisible=false)#
            else
                ax = Axis(fig[i-1, j]; xticksvisible=false, xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false)
            end

            xi = data[varnames[j]]
            yi = data[varnames[i]]

            r = round(cor(xi, yi), digits=2)
            scatter!(ax, xi, yi; kwargs...)
            scatter!(ax, 0.0, 0.0; color=(:white, 0.0), label="r = $r")
            axislegend(ax, framevisible=false)
        end
    end

    return
end
