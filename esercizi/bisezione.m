function x = bisezione(f, a, b, tol)
% x = bisezione(f, a, b, tol)
% Metodo per la ricerca degli zeri di una funzione tramite metodo di bisezione
% Input:
%   f - function che implementa la funzione f
%   a, b - estremi dell'intervallo
%   tol - tolleranza di arresto (default 1e-12)
% Output:
%   x - approssimazione della radice

if nargin < 3, error('Numero di paramametri insufficienti'), end
if nargin < 4, tol = 1e-12; end

if a >= b, error('Estremi dell'' intervallo non validi'), end
if tol <= 0, error('Tolleranza non valida'), end

fa = feval(f, a);
fb = feval(f, b);
if fa == 0, x = a; return, end
if fb == 0, x = b; return, end
if fa * fb > 0, error('Non esiste radice della funzione nell''intervallo inserito'), end


maxit = ceil(log2(b-a) - log2(tol));

for i = 1:maxit
    x = (a+b)/2;
    fx = feval(f, x);
    def = abs(fb-fa)/(b-a);
    if abs(fx) <= tol * def 
        break;
    elseif fa * fx < 0
        b = x;
        fb = fx;
    else
        a = x;
        fa = fx;
    end
end
end