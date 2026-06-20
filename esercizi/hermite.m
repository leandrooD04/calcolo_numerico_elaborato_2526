function [y, yi] = hermite(xi, fi, fi1, x)
% [y, yi] = hermite(xi, fi, fi1, x)
% Metodo di interpolazione polinoimiale tramite Hermite
% Input:
%   xi - vettore delle ascisse di interpolazione
%   fi - vettore dei valori della funzione nel punti xi
%   fi1 - vettore dei valori della derivata della funzione nei punti xi
%   x - vettore di punti di valutazione del polinomio
% Output:
%   y - valutazione del polinomio interpolante
%   yi - valutazioni della derivata prima
if nargin < 4, error('Parametri insufficienti'), end
if length(xi) ~= length(fi) || length(xi) ~= length(fi1), error('Vettori di lunghezze diverse'), end
if length(unique(xi)) ~= length(xi), error('Sono presenti ascisse non distinte'), end

n = length(xi);
m = 2*n;
z = zeros(m, 1);
Q = zeros(m, 1);

for i = 1:n
    z(2*i-1) = xi(i);
    z(2*i) = xi(i);
    Q(2*i-1) = fi(i);
    Q(2*i) = fi(i);
end

for i = m:-1:2
    if mod(i, 2) == 0
        Q(i) = fi1(i/2);
    else
        Q(i) = (Q(i)-Q(i-1)) / (z(i)-z(i-1));
    end
end

for j = 2 : m-1
    for i = m : -1 : j+1
        Q(i) = (Q(i) - Q(i-1)) / (z(i) - z(i-j));
    end
end

y = Q(m) * ones(size(x));
yi = zeros(size(x));

for i = m-1 : -1 : 1
    yi = yi .* (x - z(i)) + y;
    y  = y .* (x - z(i)) + Q(i);
end

end
