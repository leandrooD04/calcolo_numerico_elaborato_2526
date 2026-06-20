function [x, iter] = nonLinearNewton(f, x0, tol, maxit)
% x = nonLinearNewton(f, x0, tol, maxit)
% Metodo di risoluzione di sistemi non lineari con metodo di Newton
% Input:
%   f - function che implementa la funzione f
%   x0 - approssimazione iniziale
%   tol - tolleranza (default 1e-12)
%   maxit - numero di iterazioni massime (default 1000)
% Output:
%   x - approssimazione della soluzione del sistema
%   iter - numero di iterazioni eseguite
if nargin < 3, error('Numero parametri insufficiente'), end
if nargin < 5, maxit = 1000; end
if nargin < 4, tol = 1e-12; end
if maxit < 1, error('Numero di iterazioni massime non valido'), end
if tol <= 0, error('Tolleranza non valida'), end

x = x0;
flag = 1;
iter = 0;

for i = 1:maxit
    iter = iter + 1;
    x0 = x;
    [Fx, Jx] = feval(f, x);

    delta = LUPivotingSolve(Jx, -Fx);
    x = x0 + delta;
    if norm(x-x0) <= tol * (1+norm(x))
        flag = 0;
        break;
    end
end
if flag, warning('Accuratezza non raggiunta'), end
end
