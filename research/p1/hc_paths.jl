# Оценка трудоёмкости: число путей полиэдральной гомотопии (смешанный объём) для выбранных типов.
using HomotopyContinuation, JSON
const SYS = JSON.parsefile(joinpath(@__DIR__, "data", "hc", "systems.json"))

function build_system(d; plain=true)
    names = plain ? d["variables_plain"] : d["variables"]
    eqs = plain ? d["eqs_plain"] : d["eqs"]
    vars = [Variable(Symbol(n)) for n in names]
    env = Dict(Symbol(n) => v for (n, v) in zip(names, vars))
    sub(e) = e isa Symbol ? get(env, e, e) : (e isa Expr ? Expr(e.head, map(sub, e.args)...) : e)
    exprs = [eval(sub(Meta.parse(e))) for e in eqs]
    System(exprs; variables=vars)
end

ids = isempty(ARGS) ? [d["id"] for d in SYS[1:3]] : ARGS
for d in SYS
    d["id"] in ids || continue
    F = build_system(d)
    n = nvariables(F); m = length(expressions(F))
    t = @elapsed np = paths_to_track(F; start_system=:polyhedral)
    println(d["id"], " (", d["group"], "): vars ", n, ", eqs ", m, ", polyhedral paths ", np, "  [", round(t, digits=1), " s]")
    flush(stdout)
end
