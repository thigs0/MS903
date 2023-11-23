include("function_ms903.jl")
using DataFrames
using CSV
using Plots

df_treino = CSV.File("files\dado_testa", header=false, delim=" ")
df_teste = CSV.File("files\dado_treino", header=false, delim=" ")
df_treino = DataFrame(df_treino)
df_teste = DataFrame(df_teste)

global u = df_treino[:,2:3]
global v = df_treino[:,4]

scatter(u[:,1], u[:,2],v, group=v)
savefig("treino.png")

function φ(t::Float64, flag::Int=0)
	a = 0.5*π;
	f = a*atan(t)
	g = a/(1+t*t)
	if flag == 0
		return f
	elseif flag == 1
		return (f, g)
	else
		return g
	end
end

function f(x::Vector, ut)
	a1 = x[1]*ut[1] + x[2]*ut[2] + x[5]
	b1,c1 = φ(a1, 1)
	a2 = x[3]*ut[1]+x[4]*ut[2]+x[6];
	b2,c2 = φ(a2, 1)
	a = x[7]*b1 + x[8]*b2 + x[9];
	b,c = φ(a, 1);
	f=b;
	g = zeros(9)
	g[1] = x[7]*c1*ut[1];
	g[2] = x[7]*c1*ut[2];
	g[3] = x[8]*c2*ut[1];
	g[4] = x[8]*c2*ut[2];
	g[5] = x[7]*c1;
	g[6] = x[8]*c2;
	g[7] = b1;
	g[8] = b2;
	g[9] = 1
	g = c*g
	return f,g
end

function F(x::Vector, flag::Int=0)
	m = length(x);
  fv=0; g=zeros(m);
	for j in 1:m
		a,b = f(x, u[j,:]);
		aux = a-v[j];
		fv += aux*aux;
		g += aux*b;
	end
  fv /= 2*m
  fv += lambda*norm(x)^2
  g = g./m
	if flag == 0
		return fv;
	elseif flag == 1
    return (fv,g);
	else
		return g;
	end
end

rand.seed(9)
x0 = rand(9)
println(x)
l = Vector([0, 1e-3, 1e-2, 1e-1,1, 10, 10^2])
for i in l
    global lambda = i
    y = minimizador_lucio_grad(x0, F, 100)
end

scatter!(length(l), y)
savefig("Dados_lambdas.png")


global u = df_teste[:,2:3]
global v = df_teste[:,4]

global lambda = l[argminimum(y)]

scatter(u[:,1], u[:,2],v, group=v)
savefig("teste.png")
