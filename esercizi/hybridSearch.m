function x = hybridSearch(f, f1, a, b, delta, tol, maxit)
% x = hybridSearch(f, f1, a, b, delta, tol, maxit)
% Metodo di ricerca degli zeri di una funzione tramite un algoritmo ibrido
% tra metodo di bisezione e metodo di Newton
% Input:
%   f - function che implementa la funzione f
%   f1 - function che implementa la derivata della funzione f
%   a, b - estremi dell'intervallo
%   delta - ampiezza dell'intervallo minima nel metodo di bisezione
%   (default 1e-1)
%   tol - tolleranza di arresto (default 1e-12)
%   maxit - numero massimo di iterazioni (default 1000)

if nargin < 4, error('Parametri insufficienti'), end
if nargin < 7, maxit = 1000, end
if nargin < 6, tol = 1e-12, end
if nargin < 5, delta = 1e-1, end
if maxit < 1, error("Numero massimo di iterazioni non valido"), end
if tol <= 0, error("Tolleranza non valida"), end
if delta <= 0, error("Delta non valido"), end
if a >= b, error("Intervallo non valido"), end

fa = feval(f, a);
fb = feval(f, b);
if fa == 0, x = a; return, end
if fb == 0, x = b; return, end
if fa * fb > 0, error('Non esiste radice della funzione nell''intervallo inserito'), end

while b-a > delta
    x = (a+b)/2;
    fx = feval(f, x);
    if fx == 0, break;
    elseif fa * fx < 0
        b = x;
    else
        a = x;
        fa = fx;
    end
end

x = (a+b)/2;
flag = 1;

for i = 1:maxit
    x0 = x;
    fx = feval(f, x);
    f1x = feval(f1, x);
    if fx == 0
        flag = 0;
        break;
    end
    if f1x == 0, error('Derivata nulla'); end
    x = x0 - (fx/f1x);
    if abs(x-x0) <= (tol * (1 + abs(x)))
        flag = 0;
        break;
    end
end
if flag, warning('Accuratezza non raggiunta'), end
end


