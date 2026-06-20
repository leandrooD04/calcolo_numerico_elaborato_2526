function x = chebyshev(a, b, n)
% x = chebyshev(a, b, n)
% Metodo di calcolo delle ascisse di Chebyshev per un polinomio
% interpolante di grado n
% Input:
%   a, b - estremi dell'intervallo di interpolazione
%   n - grado del polinomio interpolante (vengono generati n+1 nodi)
% Output:
%   x - vettore dei nodi di Chebyshev sull'intervallo [a, b]
if nargin < 3, error('Numero parametri insufficiente'), end
if a >= b, error('Intervallo non valido'), end
if n ~= fix(n) || n < 0
    error('Il grado n deve essere un intero non negativo.');
end

x = cos((2 * (n : -1 : 0) + 1) * (pi / (2 * (n + 1))));
x = (a + b) / 2 + x * (b - a) / 2;

end