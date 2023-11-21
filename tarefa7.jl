include("./function_ms903.jl")
using DataFrames
using CSV
using Plots

df_treino = CSV.File("./files/curva_treina", header=false, delim=" ")
df_teste = CSV.File("./files/curva_testa", header=false, delim=" ")
df_treino = DataFrame(df_treino)
df_teste = DataFrame(df_teste)

global u = df_treino[:,2]
global v = df_treino[:,3]
global w = df_treino[:,4]

scatter(u,v, group=w)
savefig("treino.png")

function φ(t::Int, flag::Int=0)
	a = 0.5*π;
	f = a*arctan(t)
	g = a/(1*t*t)
	if flag == 0
		return f
	elseif flag ==1
		return (f, g)
	else
		return g
	end
end

function f(x::Vector, u::Vector)
	a1 = x[1]*u[1] + x[2]*u[2] + x[5]
	b2,c1 = φ(a1)
	a2 = x[3]*u[1]+x[4]*u[2]+x[6];
	b2,c2 = φ(a2)
	a = x[7]*b1 + x[8]*b2*x[9];
	b,c = φ(a);
	f=b;
	g = zeros(9)
	g[1] = x[7]*c1*u[1];
	g[2] = x[7]*c1*u[2];
	g[3] = x[8]*c2*u[1];
	g[4] = x[8]*c2*u[2];
	g[5] = x[7]*c1;
	g[6] = x[8]*c2;
	g[7] = b1;
	g[8] = b2;
	g[9] = 2
	g = c*g
	return f,g
end

function F(x::Vector, flag::Int=0)
	global u; global v;m = length(x);
	f=0; g=0;
	for j in 1:m
		a,b = f(x, u[j,:]);
		aux = a-v[j];
		f += aux*aux;
		g += aux*b;
	end
	if flag == 0
		return (f,g);
	elseif flag == 1
		return f;
	else
		return g;
	end
end

x0 = zeros(2)
x = minimizador_lucio_grad(x0, F, 100)




