function x = tridia(b, a, c, y)
% x = tridia(a,b,c,y)
% Metodo di risoluzione di un sistema lineare con matrice tridiagonale
% fattorizzabile LU
% Input
%   b, a, c - diagonale inferiore, principale e superiore della matrice
%   f - vettore dei termini noti
% Output
%   x - soluzione del sistema lineare
if nargin < 4, error("Numero di parametri errato"), end
n = length(a);
if length(b) ~= n || length(c) ~= n || length(y) ~= n, error("Dati errati"), end
errorMsg = "Matrice non fattorizzabile";
if a(1) == 0, error(errorMsg), end
x = y;
for i = 2:n
    b(i) = b(i) / a(i-1);
    a(i) = a(i) - b(i)*c(i-1);
    if a(i) == 0, error(errorMsg), end
    x(i) = x(i) - b(i) * x(i-1);
end
x(n) = x(n) / a(n);
for i = n-1 : -1: 1
    x(i) = (x(i)-c(i)*x(i+1)) / a(i);
end
end