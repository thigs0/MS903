include("./function.jl")
using DataFrames
using CSV
using Plots

df = CSV.File("./files/circulo.csv", header=false, delim=" ")
df = DataFrame(df)
u = df[:,2]
v = df[:,3]

function circ(p::Vector, R::Int)
  #p é um vetor com duas componentes(x,y)
  #R é o raio do circulo
  

  
end
scatter(u,v)
savefig("teste.png")
