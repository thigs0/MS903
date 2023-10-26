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

function w_circ(x1::Vector, x2::Vector, x3::Vector)
 # Com os parâmetros do circulo, retorna o vetor que prevemos 
  c = (x3[1].-x1).^2 + (x3[2].-x2).^2 
  c = c .- x3[3]^2
  c = sign.(c)
end

x = Vector([0.0, 0, 1]) # cada coordenada representa eixo x, eixo y e raio
x = minimizador_lucio(x, circ, 100)
circulo(x) = x[3]*sin.(range(0, 2π, 100)).+x[1], x[3]*cos.(range(0, 2π, 100)).+x[2]
plot(circulo(x))
w2 = w_circ(u,v,x);
println("A matrix de confusão dos dados de treino é:")
display(confusion_matrix(w,w2))

scatter!(u,v, group=w)
savefig("dados_treino.png")

global u = df_teste[:,2]
global v = df_teste[:,3]
global w = df_teste[:,4]
w2 = w_circ(u,v,x);

plot(circulo(x))
scatter!(u,v, group=w)
println("A matrix de confusão dos dados de teste é:")
display(confusion_matrix(w,w2))
savefig("dados_teste.png")
