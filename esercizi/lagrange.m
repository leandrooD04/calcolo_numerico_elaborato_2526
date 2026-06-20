function y = lagrange(xi, fi, x)
% x = lagrange(xi, fi, x)
% Metodo di interpolazione polinomiale in forma di Lagrange
% Input:
%   xi - vettore delle ascisse di interpolazione
%   fi - vettore dei valori della funzione nel punti xi
%   x - vettore di punti di valutazione del polinomio
% Output:
%   y - valutazione del polinomio interpolante
if nargin < 3, error('Numero parametri insufficiente'), end
if length(xi) ~= length(fi), error('Vettori di lunghezze diverse'), end
if length(unique(xi)) ~= length(xi), error('Sono presenti ascisse non distinte'), end

n = length(xi);
y = zeros(size(x));

for i = 1:n
    Li = ones(size(x));

    for j = 1:n
        if i ~= j
            Li = Li .* (x - xi(j)) / (xi(i) - xi(j));
        end
    end
    y = y + fi(i) * Li;
end
end