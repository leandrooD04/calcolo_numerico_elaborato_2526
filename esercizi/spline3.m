function y = spline3(xi, fi, x, type, f1a, f1b)
% y = spline3(xi, fi, x, type, f1a, f1b)
% Spline cubica interpolante via metodo dei momenti (sistema tridiagonale)
% Input:
%   xi - ascisse di interpolazione
%   fi - valori della funzione nei punti xi
%   x  - punti dove valutare la spline
%   type - 0 naturale, 1 not-a-knot (default), 2 completa
%   f1a, f1b - derivate prime agli estremi (richieste solo se type==2)
% Output:
%   y - valutazioni della spline (stessa forma di x)

if nargin < 3, error('Parametri insufficienti'), end
if nargin < 4, type = 1; end
if type == 2 && nargin < 6, error('Per la spline completa servono f1a e f1b'), end

xi = xi(:);
fi = fi(:);
n = length(xi) - 1;
if n + 1 ~= length(fi), error('Lunghezze dei vettori in ingresso non valide'), end
if n < 2, error('Servono almeno 3 nodi'), end
if type == 1 && n < 3, error('La spline not-a-knot richiede almeno 4 nodi'), end

h     = diff(xi);
delta = diff(fi) ./ h;

phi = zeros(n+1,1);  ksi = zeros(n+1,1);  dd = zeros(n+1,1);
for k = 2:n
    hl = h(k-1);  hr = h(k);
    phi(k) = hl / (hl + hr);
    ksi(k) = hr / (hl + hr);
    dd(k)  = 6*(delta(k) - delta(k-1)) / (xi(k+1) - xi(k-1));
end

m = zeros(n+1, 1);

if type == 0 || type == 2
    b = zeros(n+1,1);  a = zeros(n+1,1);  c = zeros(n+1,1);  rhs = zeros(n+1,1);
    a(2:n) = 2;
    b(2:n) = phi(2:n);
    c(2:n) = ksi(2:n);
    rhs(2:n) = dd(2:n);
    if type == 0    % NATURALE
        a(1) = 1;  c(1) = 0;  rhs(1) = 0;
        a(n+1) = 1;  b(n+1) = 0;  rhs(n+1) = 0;
    else            % COMPLETA
        a(1) = 2;  c(1) = 1;  rhs(1) = 6*(delta(1) - f1a) / h(1);
        a(n+1) = 2;  b(n+1) = 1;  rhs(n+1) = 6*(f1b - delta(n)) / h(n);
    end
    m = tridia(b, a, c, rhs);

else                % NOT-A-KNOT
    b = zeros(n-1,1);  a = zeros(n-1,1);  c = zeros(n-1,1);  rhs = zeros(n-1,1);
    for j = 1:n-1
        k = j + 1;
        b(j) = phi(k);  a(j) = 2;  c(j) = ksi(k);  rhs(j) = dd(k);
    end
    a(1) = 2 - phi(2);  c(1) = ksi(2) - phi(2);  rhs(1) = ksi(2) * dd(2);
    b(n-1) = phi(n) - ksi(n);  a(n-1) = 2 - ksi(n);  rhs(n-1) = phi(n) * dd(n);

    mInt = tridia(b, a, c, rhs);
    m(2:n) = mInt;
    m(1)   = dd(2) - m(2)   - m(3);
    m(n+1) = dd(n) - m(n-1) - m(n);
end

y = zeros(size(x));
for j = 1:length(x)
    if x(j) <= xi(1)
        i = 1;
    elseif x(j) >= xi(n+1)
        i = n;
    else
        i = find(xi <= x(j), 1, 'last');
        if i > n, i = n; end
    end
    d1 = xi(i+1) - x(j);
    d0 = x(j) - xi(i);
    y(j) = (m(i)*d1^3 + m(i+1)*d0^3)/(6*h(i)) ...
         + (fi(i)/h(i)   - h(i)*m(i)/6)   * d1 ...
         + (fi(i+1)/h(i) - h(i)*m(i+1)/6) * d0;
end
end