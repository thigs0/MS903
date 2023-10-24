include("./function_ms903.jl")

using DataFrames
using CSV
using Plots

df_treino = CSV.File("./files/circulo_treina", header=false, delim=" ")
df_teste = CSV.File("./files/circulo_testa", header=false, delim=" ")
df_treino = DataFrame(df_treino)
df_teste = DataFrame(df_teste)

global u = df_treino[:,2]
global v = df_treino[:,3]
global w = df_treino[:,4]

function circ(x::Vector)
  #p é um vetor com duas componentes(x,y)
 # é um vetor com o centro da esfera
  #R é o raio do circulo
  # A função diz se temos um ponto dentro ou fora da circunferencia
  
  c = (x[1].-u).^2 + (x[2].-v).^2;
  c = c .-x[3]^2
  aux1 = atan.(c);
  return sum((aux1 - (pi/2)*w).^2)/length(u);
end

x = Vector([0.0, 0, 1])
x = minimizador_lucio(x, circ, 100)
tmin = 0
tmax = 2π
tvec = range(tmin, tmax, length = 100)

plot(x[3]*sin.(tvec).+x[1], x[3]*cos.(tvec).+x[2])

scatter!(u,v, group=w)
savefig("teste.png")
