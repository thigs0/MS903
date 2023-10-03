include("./function_ms903.jl")

#abre os dados e capta o u, v
using CSV
using Plots
using DataFrames
using LinearAlgebra

csv = CSV.File("./files/tarefa3")
csv = DataFrame(csv)

u = csv[:, 2];
v = csv[:, 3];
w = csv[:, 4];

x = [7; 4];
r = grad_est(x, u, v, w, 41)
s = scatter(u,v, group=w)
f(x) = r[2]*x+r[1]
plot!(f, -10, 10)
savefig(s, "save.png")

x = [7; 4];
r = grad_est_p(x, u, v, w, 100, 41)
s = scatter(u,v, group=w)
f(x) = r[2]*x+r[1]
plot!(f, -10, 10)
savefig(s, "save.png")
