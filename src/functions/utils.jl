function includedir(dir)
    include.(filter(contains(r".jl$"), readdir(dir; join = true)))
end


function readdirext(dir, ext; join = false)
    dc = readdir(dir; join = join)
    out = String[]
    N = length(ext) - 1
    for d in dc
        if d[(end-N):end] == ext
            push!(out, d)
        end
    end
    return out
end
