function mesh2obj(m, outfile)
    open(outfile, write=true) do io
        for v in coordinates(m)
            write(io, "v $(join(v, " "))\n")
        end
        for f in faces(m)
            write(io, "f $(join(f, " "))\n")
        end
    end
end
