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
            xticksvisible = false,
            yticksvisible = false,
            xlabelpadding = 3,
            ylabelpadding = 3
        ),
        Legend = (
            framevisible = false,
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
    θ = LinRange(0, π, n)
    φ = LinRange(-π, π, 2n)
    xe = [cos(φ)*sin(θ) for θ in θ, φ in φ]
    ye = [sin(φ)*sin(θ) for θ in θ, φ in φ]
    ze = [cos(θ) for θ in θ, φ in φ]
    surface!(ax, R*xe .+ pos[1], R*ye .+ pos[2], R*ze .+ pos[3];
        color=load(joinpath(@__DIR__, "../../assets/lroc_color_poles_2k.tif")), shading=false)
end

function axisoff!(ax)
    hidedecorations!(ax)  # hides ticks, grid and lables
    hidespines!(ax)  # hide the frame
end

# Plot axes of frame F in world (image) frame W
function plotframe!(ax, pos_W, R_WF, length=1)
    xF_W = R_WF[:, 1]
    yF_W = R_WF[:, 2]
    zF_W = R_WF[:, 3]
    arrows!(ax, [pos_W[1]], [pos_W[2]], [pos_W[3]], [xF_W[1]], [xF_W[2]], [xF_W[3]]; arrowsize=0.1*length, lengthscale=length, color=:red)
    arrows!(ax, [pos_W[1]], [pos_W[2]], [pos_W[3]], [yF_W[1]], [yF_W[2]], [yF_W[3]]; arrowsize=0.1*length, lengthscale=length, color=:blue)
    arrows!(ax, [pos_W[1]], [pos_W[2]], [pos_W[3]], [zF_W[1]], [zF_W[2]], [zF_W[3]]; arrowsize=0.1*length, lengthscale=length, color=:green)
end

# Transform 3D model
function transformModel(m, pos_I, R_IB=I, scale=1.0)
    c, f = coordinates(m), faces(m)
    mt = [Point3f(pos_I + scale*R_IB*c[k]) for k in eachindex(c)]
    return GeometryBasics.Mesh(mt, f)
end
