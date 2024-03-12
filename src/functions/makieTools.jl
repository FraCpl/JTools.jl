function theme_fra(spinevis=true)
    bg = RGBf(40/255, 44/255, 52/255)
    fg = RGBf(0.7,0.7,0.7)
    fga = RGBf(0.3,0.3,0.3)
    Theme(
        backgroundcolor = bg,
        textcolor = fg,
        linecolor = fga,
        Axis = (
            backgroundcolor = :transparent,
            xgridcolor = (:white, 0.09),
            ygridcolor = (:white, 0.09),
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
            xminorgridcolor = (:white, 0.02),
            yminorgridcolor = (:white, 0.02),
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
            xgridcolor = (:white, 0.09),
            ygridcolor = (:white, 0.09),
            zgridcolor = (:white, 0.09),
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
#ax = Axis3(f[1, 1]; aspect = (1, 1, 1))
#plotMoon!(ax)
function plotMoon!(ax, pos=zeros(3), R=1.0; n=256)
    θ = LinRange(0, π, n)       # Colatitude
    φ = LinRange(-π, π, 2n)     # Longitude
    xe = [cos(φ)*sin(θ) for θ in θ, φ in φ]
    ye = [sin(φ)*sin(θ) for θ in θ, φ in φ]
    ze = [cos(θ) for θ in θ, φ in φ]
    surface!(ax, R*xe .+ pos[1], R*ye .+ pos[2], R*ze .+ pos[3];
        color=load(joinpath(@__DIR__, "../../assets/lroc_color_poles_2k.tif")), shading=NoShading)
end

function axisoff!(ax)
    hidedecorations!(ax)  # hides ticks, grid and lables
    hidespines!(ax)  # hide the frame
end

# Plot axes of frame F in world (image) frame W
# CAUTION: This seems to make MAKIE crash
function plotframe!(ax, pos_W, R_WF, length=1)
    xF_W = R_WF[:, 1]*length
    yF_W = R_WF[:, 2]*length
    zF_W = R_WF[:, 3]*length
    # CAUTION: This seems to make MAKIE crash
    #arrows!(ax, [pos_W[1]], [pos_W[2]], [pos_W[3]], [xF_W[1]], [xF_W[2]], [xF_W[3]]; arrowsize=0.1*length, lengthscale=length, color=:red)
    #arrows!(ax, [pos_W[1]], [pos_W[2]], [pos_W[3]], [yF_W[1]], [yF_W[2]], [yF_W[3]]; arrowsize=0.1*length, lengthscale=length, color=:blue)
    #arrows!(ax, [pos_W[1]], [pos_W[2]], [pos_W[3]], [zF_W[1]], [zF_W[2]], [zF_W[3]]; arrowsize=0.1*length, lengthscale=length, color=:green)
    lines!(ax, pos_W[1] .+ [0.0, xF_W[1]], pos_W[2] .+ [0.0, xF_W[2]], pos_W[3] .+ [0.0, xF_W[3]]; color=:red)
    lines!(ax, pos_W[1] .+ [0.0, yF_W[1]], pos_W[2] .+ [0.0, yF_W[2]], pos_W[3] .+ [0.0, yF_W[3]]; color=:blue)
    lines!(ax, pos_W[1] .+ [0.0, zF_W[1]], pos_W[2] .+ [0.0, zF_W[2]], pos_W[3] .+ [0.0, zF_W[3]]; color=:green)
end

# Transform 3D model
function transformModel(m, pos_I, R_IB=I, scale=1.0)
    c, f = coordinates(m), faces(m)
    mt = [Point3f(pos_I + scale*R_IB*c[k]) for k in eachindex(c)]
    return GeometryBasics.Mesh(mt, f)
end

function sensorFovModel(; pos_I=zeros(3), R_IS=I, FOV=(40π/180, 40π/180), length=1)
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

function cuboidModel(lx, ly, lz; pos_I=zeros(3), R_IB=I)
    v = [Point3f(pos_I + R_IB*[+lx/2; -ly/2; -lz/2]),
         Point3f(pos_I + R_IB*[+lx/2; +ly/2; -lz/2]),
         Point3f(pos_I + R_IB*[+lx/2; +ly/2; +lz/2]),
         Point3f(pos_I + R_IB*[+lx/2; -ly/2; +lz/2]),
         Point3f(pos_I + R_IB*[-lx/2; -ly/2; -lz/2]),
         Point3f(pos_I + R_IB*[-lx/2; +ly/2; -lz/2]),
         Point3f(pos_I + R_IB*[-lx/2; +ly/2; +lz/2]),
         Point3f(pos_I + R_IB*[-lx/2; -ly/2; +lz/2]),
    ]
    f = [QuadFace([1, 2, 3, 4]), QuadFace([2, 3, 7, 6]), QuadFace([7, 8, 5, 6]), QuadFace([1, 4, 8, 5]), QuadFace([4, 3, 7, 8]), QuadFace([1, 2, 6, 5])]
    return GeometryBasics.Mesh(v, f)
end

function nicholsgrid(f; center::Int=-1)
    ct = 360*center + 180
    ax = Axis(f; limits=(ct - 180, ct + 180, -40, 40), xgridvisible=false, ygridvisible=false,
        xlabel="Open-Loop Phase [deg]", ylabel="Open-Loop Gain [dB]", title="Nichols Chart")
    phase = range(1, 359, 1000)*π/180;

    for mag in [6; 3; 2; 1; 0.5; 0; -0.5; -1; -2; -3; -6; -9; -12; -15; -20; -25; -30; -40]
        H = db2mag(mag).*exp.(im*phase)
        HH = H./(1 .- H)
        ph = atan.(imag.(HH), real.(HH))
        unwrap!(ph)
        mg = mag2db.(abs.(HH))

        lines!(ax, ph*180/pi .+ (ct - 180), mg, linewidth=0.5, color=:grey, linestyle=:dash)
    end

    mag = -40:0.05:6
    for phs = [359; 355; 350:-10:10; 5; 1]
        H = db2mag.(mag).*exp(im*phs*pi/180)
        HH = H./(1 .- H)
        ph = atan.(imag.(HH), real.(HH))
        mg = mag2db.(abs.(HH))

        lines!(ax, ph*180/pi .+ 360.0.*(ph .≤ 0.0) .+ (ct - 180), mg, linewidth=0.5, color=:grey, linestyle=:dash)
    end
    scatter!(ax, ct, 0, marker=:cross, color=:grey, markersize=8)
    return ax
end
