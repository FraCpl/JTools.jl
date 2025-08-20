function mesh2obj(m, outfile, flipNormal=false)
    open(outfile, write=true) do io
        for v in coordinates(m)
            write(io, "v $(join(v, " "))\n")
        end
        if flipNormal
            flpidx = [1, 3, 2]
            for f in faces(m)
                write(io, "f $(join(f[flpidx], " "))\n")
            end
        else
            for f in faces(m)
                write(io, "f $(join(f, " "))\n")
            end
        end
    end
end

# Pos is a matrix of 3d positions
function grid2mesh(pos)
    Ny, Nx = size(pos)
    f = Vector{TriangleFace}(undef, 2(Nx - 1)*(Ny - 1))
    q = 1
    for j in 1:Nx-1, k in Ny*(j - 1)+1:Ny*j-1
        f[q] = TriangleFace(k, k + 1, k + Ny + 1)
        f[q+1] = TriangleFace(k, k + Ny + 1, k + Ny)
        q += 2
    end
    return GeometryBasics.Mesh(Point3f.(pos[:]), f)
end

struct ObjModel
    faces::Vector{Vector{Int}}
    vertices::Vector{Vector{Float64}}
    edges::Vector{Vector{Int}}
    normals::Vector{Vector{Float64}}
end

function ObjModel(faces, vertices; computeNormals=false, computeEdges=false)
    # Compute model edges
    e = Vector{Vector{Int}}(undef, 0)
    if computeEdges; e = compEdges(faces); end

    # Compute model normals
    n = Vector{Vector{Float64}}(undef, 0)
    if computeNormals; n = compNormals(faces, vertices); end

    return ObjModel(faces, vertices, e, n)
end

function writeObj(obj::ObjModel, outfile::String)
    open(outfile, write=true) do io
        for v in obj.vertices
            write(io, "v $(join(v, " "))\n")
        end
        for f in obj.faces
            write(io, "f $(join(f, " "))\n")
        end
    end
end

function readObj(objfile::String; computeNormals=false, computeEdges=false)
    v = Vector{Vector{Float64}}(undef, 0)
    f = Vector{Vector{Int}}(undef, 0)
    open(objfile) do fl
        while !eof(fl)
            s = readline(fl)
            if length(s) > 2
                if s[1:2] == "v "
                    push!(v, parse.(Float64, split(s," "; keepempty=false)[2:4]))
                end
                if s[1:2] == "f "
                    fk = Vector{Int}(undef, 0)
                    for fc in split(s," "; keepempty=false)[2:end]
                        push!(fk, parse(Int, split(fc, "/"; keepempty = false)[1]))
                    end
                    push!(f, fk)
                end
            end
        end
    end

   # Compute model edges
   e = Vector{Vector{Int}}(undef, 0)
   if computeEdges; e = compEdges(f); end

   # Compute model normals
   n = Vector{Vector{Float64}}(undef, 0)
   if computeNormals; n = compNormals(f, v); end

    return ObjModel(f, v, e, n)
end

function compEdges(f)
    e = Vector{Vector{Int}}(undef, 0)
    for face in f
        fi = sort(face)
        e1 = [fi[1], fi[2]]; e2 = [fi[2], fi[3]]; e3 = [fi[1], fi[3]]
        if !any(Ref(e1) .== e); push!(e, e1); end
        if !any(Ref(e2) .== e); push!(e, e2); end
        if !any(Ref(e3) .== e); push!(e, e3); end
    end
    return e
end

function compNormals(f, v)
    E1 = [v[face[2]] - v[face[1]] for face in f]
    E2 = [v[face[3]] - v[face[2]] for face in f]
    return normalize.(E1 .× E2)
end
