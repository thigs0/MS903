# MS903

Conteúdo para fixar o aprendizado de máquina do vsito matemático com o como a função irá aproximar a sua saída com o que é esperado.

### Problema inverso
Dado uma função f(x) = y com as saídas, muitas vezes temos a função e a saída, a fim de descobrir o ponto que originou. Já quando temos o x e a função, queremos saber calacular a saída de forma otimizada.

Em nosso problema, temos o ponto x e o ponto y. Assim queremos construir uma função geral que aproxime os dados sem saber como a função se comporta. 
Não confunda isso com a otimização pura de parâmetros de uma função como uma regrassão linear, visto que nesse caso já temos a função. Imagine que não temos a função.

A princípio criaremos uma função simples apenas de adição e multiplicação crescente até conseguirmos um bom resultado.


Criaremos as funçãos que usaremos a frente

´´´

function minimizador(x::Vector{Float64}, f::Function, M::Int, α=10e-4, σ=0.5, ε=10e-9) 
    #=Dado um vetor inicial, uma função que é aplicada em x 
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

´´´
