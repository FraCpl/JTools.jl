function includedir(dir)
    include.(filter(contains(r".jl$"), readdir(dir; join=true)))
end
