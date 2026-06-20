function [Ikn, err] = ncotescomp(fx, a, b, k, n)
% [Ikr, err] = ncotescomp(fx, a, b, k, n)
% Metodo di calcolo dell'approssimazione dell'integrale di una funzione f 
% nell'intervallo [a, b] tramite la formula composita di Newton-Cotes
% Input:
%   fx - function che implementa la funzione f
%   a, b - estremi di integrazione
%   k - grado della formula di Newton-Cotes
%   n - numero di suddivisioni dell'intervallo
% Output:
%   Ikn - approssimazione dell'integrale
%   err - stima dell'errore di quadratura
if nargin < 5, error('Numero parametri insufficiente'), end
if a >= b, error('Intervallo di integraizone non valido'), end
if k <= 0 || k ~= fix(k), error('Grado della formula non valido'), end
if n <= 0 || mod(n, k) ~= 0 || mod(n/k, 2) ~= 0, error('Numero di suddivisioni dell''intervallo non valido'), end

h = (b - a) / n;
x = linspace(a, b, n+1);
f = feval(fx, x);
f = f(:);
ci = ncotes(k);
ci = ci(:);

Ikn = 0;
for i = 1:n/k
    idx = (i-1)*k+1:i*k+1;
    Ikn = Ikn + h * sum(ci .* f(idx));
end

Ikn2 = 0;
for i = 1:n/(2*k)
    idx = (i-1)*2*k+1:2:i*2*k+1;
    Ikn2 = Ikn2 + (2*h) * sum(ci .* f(idx));
end

if mod(k, 2) == 0
    p = k+2;
else
    p = k+1;
end

err = (Ikn - Ikn2) / (2^p - 1);
end