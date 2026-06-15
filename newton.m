function x = newton(f, f1, x0, tol, maxit, molt)
% function x = newton(f, f1, x0, tol, maxit, molt)
% metodo per la ricerca di zeri di una funzione tramite metodo di Newton
% modificato
% Input:
%   f - function che implementa la funzione f
%   f1 - function che implementa la derivata della funzione f
%   x0 - punto iniziale
%   tol - tolleranza di arresto (default = 1e-12)
%   maxit - numero massimo di iterazioni (default = 1000)
%   molt - molteplicita' della radice (default = 1)
% Output:
%   x - approssimazione della radice

if nargin < 3, error('Numero insufficiente di parametri'), end
if nargin < 6, molt = 1; end
if nargin < 5, maxit = 1000; end
if nargin < 4, tol = 1e-12; end

if molt < 1 || molt ~= fix(molt), error('Molteplicita'' non valida'); end
if maxit <= 1, error('Numero di interazioni massime non valido'); end
if tol <= 0, error('Tolleranza non valida'); end

x = x0; 
flag = 1;
for i=1:maxit
    x0 = x;
    fx = feval(f, x);
    f1x = feval(f1, x);
    if fx == 0
        flag = 0;
        break;
    end
    if f1x == 0, error('Derivata nulla'); end
    x = x0 - molt * (fx/f1x);
    if abs(x-x0) <= (tol * (1 + abs(x)))
        flag = 0;
        break;
    end
end
if flag, warning('Accuratezza non raggiunta'), end
end