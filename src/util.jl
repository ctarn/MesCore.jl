using Dates

match_path(path, ext="") = begin
    paths = readdir(dirname(path); join=true)
    return filter(p -> startswith(basename(p), basename(path)) && endswith(p, ext), paths)
end

read_all(reader, path, ext=""; single=false, key=first∘splitext∘basename) = begin
    paths = match_path(path, ext)
    if single && length(paths) > 1
        @warn "multiple files found: path=$(path), ext=$(ext)"
    end
    return map(paths) do p
        @info "reading from " * p
        key(p) => reader(p)
    end |> Dict
end

parse_range(::Type{T}, s::AbstractString) where T = begin
    vs = map(x -> parse(T,  x),  split(s, ':'))
    if length(vs) <= 2 return range(vs[begin], vs[end])
    elseif length(vs) == 3 return range(vs[begin], vs[end]; step=vs[begin+1])
    else return nothing
    end
end

date_mark(date=today()) = replace(string(date), "-"=>"")

open_url(url) = begin
    if Sys.isapple()
        run(`open $(url)`)
    elseif Sys.iswindows()
        run(`cmd /k start $(url) "&&" exit`)
    else
        run(`open $(url)`)
    end
end
