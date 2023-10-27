using Random
using LinearAlgebra
using Calculus

function grad_partial(x, u, v, w)
  #=
  
	=#
  ale = floor(Int64, length(u)*rand())+1 # numero aleatório de 1 até n
  aux1 = x[1] + x[2]*u[ale]-v[ale]
  aux2 = atan(aux1)+0.5*pi*w[ale];
  aux3 = aux2/(1 + aux1*aux1);
  g = zeros(2)

  g[1] = aux3
  g[2] = aux3*u[ale]
  return g
end

function grad_partial_p(x, u, v, w, p)
  #=
  
	=#
  g = zeros(2)

  for i in rand(1:length(u), p)
    aux1 = x[1] + x[2]*u[i]-v[i]
    aux2 = atan(aux1)+0.5*pi*w[i];
    aux3 = aux2/(1 + aux1*aux1);
    g[1] += aux3
    g[2] += aux3*u[i]
  end

  return g
end

function quadratica(x::Vector{Float64})
  #=
  #Calcula a função quadrática para o vetor x e retorna um float
  =#
    y = 0
    x = x.^2
    for i in 1:length(x)
        y += i*x[i]
    end
    y /= 2
    return y
end

function Rosenbrock(x::Vector{Float64})
  #=
  #Calcula a função de Rosenbrock para o vetor x e retorna um float
  =#
    y=0
    for i in 1:Int((length(x)/2))
        y += 10*(x[2*i]-x[2*i-1]^2)^2 + (x[2*i-1]-1)^2
    end
    return y
end

function fun(x::Vector{Float64}, grad=0)
	#=grad(Int) -1 retorna o gradiente
				0 retorna a função 
				1 retorna a função e o gradiente
	=#
	aux1 = x[1] .+ x[2].*u-v
	aux2 = atan.(aux1)+0.5*pi*w;
  
	f = 0.5*sum(aux1.*aux2)
  if grad == 1
	  aux3 = aux2./(1 .+ aux1.*aux1);
    g = zeros(2)
  	g[1] = sum(aux3)/n
	  g[2] = sum(aux3.*u)/n
	  return f, g
	elseif grad == 0
	    return f;
	else 
		aux3 = aux2./(1 .+ aux1.*aux1);
	    g = zeros(2)
	  	g[1] = sum(aux3)/n
		g[2] = sum(aux3.*u)/n
		return g
  end
end


function minimizador(x::Vector{Float64}, f::Function, M::Int, α=10e-4, σ=0.5, ε=10e-9)
    x = copy(x)
    x_p = copy(x)
    k=0; g= Calculus.gradient(f, x_p);
    g_p = copy(g);
    η = norm(g_p, Inf) 

    while (η >= ε) && (k < M)
        if k==0
            t=1
        else
            t = norm(x_p-x, 2)/norm(g_p-g, 2)
        end
        W = α*(g_p')*g_p;
        fx = f(x_p)
        while f(x_p-t*g_p) > fx-t*W
            t =  σ*t;
        end
        x = copy(x_p);
        x_p = x_p-t*g_p;
        g = copy(g_p)
        g_p = Calculus.gradient(f, x_p)
        η = norm(g_p, Inf);
        k = k+1
    end
    return x_p;
end

function minimizador_lucio_grad(x::Vector{Float64}, f::Function, M::Int, α=10e-4, σ=0.5, ε=10e-9)
    x = copy(x0)
    x_p = copy(x)
    k=0; g= f(x, -1)
    g_p = copy(g);
    η = norm(g_p, Inf) 

    while (η >= ε) && (k < M)
        if k==0
            t=1
        else
            t = norm(x_p-x, 2)/norm(g_p-g, 2)
        end
        w = α*(g_p')*g_p;
        fx = f(x_p)
        while f(x_p-t*g_p) > fx-t*w
            t =  σ*t;
        end
        x = copy(x_p);
        x_p = x_p-t*g_p;
        g = copy(g_p)
        g_p = f(x_p, -1)
        η = norm(g_p, Inf);
        k = k+1
    end
    return x_p;
end

function grad_est(x::Vector{Float64}, u, v, w, M::Int, t::Int=1, ε=10e-9)
    x = copy(x);
    x_p = copy(x);
    g = grad_partial(x_p, u, v, w); k=0; 
    g_p = copy(g);
    η = norm(g, Inf);

    while (η >= ε) && (k < M)

        x = copy(x_p);
        x_p = x_p-t*g_p;
        g = copy(g_p)
        g_p = grad_partial(x_p, u, v, w)
        η = norm(g_p, Inf);
        k+=1
    end
    return x_p;
end

function grad_estp(x::Vector, u, v, w ,p::Int , M, t=1, ε=10e-9)
    x = copy(x);
    x_p = copy(x);
    g = grad_partial_p(x_p, u, v, w); k=0; 
    g_p = copy(g);
    η = norm(g, Inf);

    while (η >= ε) && (k < M)

        x = copy(x_p);
        x_p = x_p-t*g_p;
        g = copy(g_p)
        g_p = grad_partial_p(x_p, u, v, w, p)
        η = norm(g_p, Inf);
        k+=1
    end
    return x_p;
end

function confusion_matrix(w1::Vector, w2::Vector)
  #= Calcula a matriz de confusão sendo que w1 são os valores exatas e w2 são os previstos pelo modelo
  =#
  M = [0.0 0.0; 0.0 0.0]; # matriz de confusão
  n = length(w1)
  for i in 1:n
    if w1[i] == w2[i] # Se acertamos na previsão
      if w1[i] == 1
        M[1, 1] += 1
      else
        M[2, 2] += 1
      end

    else # Se erramos a previsão
      if w1[i] == 1 # se o esperado era positivo
        M[1,2] += 1
      else
        M[2,1] += 1
      end
    end
  end
  return M
end
