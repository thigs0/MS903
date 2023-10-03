using Random
using Calculus
using LinearAlgebra

module Gradientes 
  function grad_partial(x::Vector, u::Vector, v::Vector, w::Vector)
    #= A função calcula o gradiente estocastico para os vetores (u, v) classificados por w
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

  function grad_partial_p(x::Vector, u::Vector, v::Vector, w::Vector, p::Int64)
    #= A função calcula o gradiente estocastico para os vetores (u, v) classificados por w
      usando p componetes
    
    =#
    g = zeros(2)
    if p > length(u)
      error("o valor de p precisa ser menor do que a quantidade de componentes de u")
    end
    for i in rand(1:length(u), p)
      aux1 = x[1] + x[2]*u[i]-v[i]
      aux2 = atan(aux1)+0.5*pi*w[i];
      aux3 = aux2/(1 + aux1*aux1);
      g[1] += aux3
      g[2] += aux3*u[i]
    end
    return g
  end
end

module funcoes
  function quadratica(x::Vector)
      #=
    #Calcula a função quadrática para o vetor x
    =#
      y = 0
      x = x.^2
      for i in 1:length(x)
          y += i*x[i]
      end
      y /= 2
      return y
  end
  function Rosenbrock(x::Vector)
    #= Calcula a função de Rosenbrock para o vetor x
    =#
      y=0
      for i in 1:Int((length(x)/2))
          y += 10*(x[2*i]-x[2*i-1]^2)^2 + (x[2*i-1]-1)^2
      end
      return y
  end
  function fun(x::Vector,u::Vector, v::Vector,w::Vector, grad::Int64=0)
    #=grad(Int) -1 retorna o gradiente
          0 retorna a função 
          1 retorna a função e o gradiente
    =#
    aux1 = x[1] .+ x[2].*u-v
    aux2 = atan.(aux1)+0.5*pi*w;
    n = length(u)
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
end

function minimizador_lucio(x0::Vector, f::Function, M::Int64, α::Int64=10e-4, σ::Int64=0.5, ε::Int64=10e-9)
    x_ant = copy(x0)
    x_pos = copy(x0)
    k=0; g_ant= Calculus.gradient(f, x_pos);
    g_pos = copy(g_ant);
    η = norm(g_pos, Inf) 
    while (η >= ε) && (k < M)
        if k==0
            t=1
        else
            t = norm(x_pos-x_ant, 2)/norm(g_pos-g_ant, 2)
        end
        w = α*(g_pos')*g_pos;
        fx = f(x_pos)
        while f(x_pos-t*g_pos) > fx-t*w
            t =  σ*t;
        end
        x_ant = copy(x_pos);
        x_pos = x_pos-t*g_pos;
        g_ant = copy(g_pos)
        g_pos = Calculus.gradient(f, x_pos)
        η = norm(g_pos, Inf);
        k = k+1
    end
    return x_pos;
end

function minimizador_lucio_grad(x0, f, M, α=10e-4, σ=0.5, ε=10e-9)
    x_ant = copy(x0)
    x_pos = copy(x0)
    k=0; g_ant= f(x_ant, -1)
    g_pos = copy(g_ant);
    η = norm(g_pos, Inf) 

    while (η >= ε) && (k < M)
        if k==0
            t=1
        else
            t = norm(x_pos-x_ant, 2)/norm(g_pos-g_ant, 2)
        end
        w = α*(g_pos')*g_pos;
        fx = f(x_pos)
        while f(x_pos-t*g_pos) > fx-t*w
            t =  σ*t;
        end
        x_ant = copy(x_pos);
        x_pos = x_pos-t*g_pos;
        g_ant = copy(g_pos)
        g_pos = f(x_pos, -1)
        η = norm(g_pos, Inf);
        k = k+1
    end
    return x_pos;
end

function grad_est(x0::Vector, u::Vector, v::Vector, w::Vector, M::Int64, t::Int64=1, ε::Int64=10e-9)
    x = copy(x0);
    x_pos = copy(x0);
    g = grad_partial(x_pos, u, v, w); k=0; 
    g_pos = copy(g);
    η = norm(g, Inf);

    while (η >= ε) && (k < M)

        x = copy(x_pos);
        x_pos = x_pos-t*g_pos;
        g = copy(g_pos)
        g_pos = grad_partial(x_pos, u, v, w)
        η = norm(g_pos, Inf);
        k+=1
    end
    return x_pos;
end

function Minimizador_grad_estp(x0::Vector, u::Vector, v::Vector, w::Vector,p::Int64 , M::Int64, t::Int64=1, ε::Int64=10e-9)
  #= A função minimiza 
  =#
    x = copy(x0);
    x_pos = copy(x0);
    g = grad_partial_p(x_pos, u, v, w); k=0; 
    g_pos = copy(g);
    η = norm(g, Inf);
    while (η >= ε) && (k < M)
        x = copy(x_pos);
        x_pos = x_pos-t*g_pos;
        g = copy(g_pos)
        g_pos = grad_partial_p(x_pos, u, v, w, p)
        η = norm(g_pos, Inf);
        k+=1
    end
    return x_pos;
end


