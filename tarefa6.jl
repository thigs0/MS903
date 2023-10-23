include("./function.jl")
using DataFrames
using CSV
using Plots

df_treino = CSV.File("./files/circulo_treina", header=false, delim=" ")
df_teste = CSV.File("./files/circulo_testa", header=false, delim=" ")
df_treino = DataFrame(df_treino)
df_teste = DataFrame(df_teste)

u = df_treino[:,2]
v = df_treino[:,3]
w = df_treino[:,4]

function circ(x::Vector,C::Vector ,R::Int)
  #p é um vetor com duas componentes(x,y)
 # é um vetor com o centro da esfera
  #R é o raio do circulo
  # A função diz se temos um ponto dentro ou fora da circunferencia
  
  c = (x[1]-C[1])^2 + (x[2]-C[2])^2 - R^2
  return atan(c+(pi/2))
end


scatter(u,v, group=w)
savefig("teste.png")
