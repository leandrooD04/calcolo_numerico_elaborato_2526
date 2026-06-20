function ci = ncotes(n)
% ci = ncotes(n)
% Metodo di calcolo dei coefficienti della formula di Newton-Cotes di
% grado n
% Input:
%   n - grado della formula
% Output:
%   ci - vettore coefficienti di grado n
if nargin < 1, error('Numero parametri insuficiente'), end
if n <= 0, error('Grado della formula non valido'), end

ci = zeros(n+1, 1);
for i = 0:n
    ind = [0:i-1, i+1:n];
    din = prod(i-ind);
    a = poly(ind);
    q = [a./(n+1:-1:1), 0];
    alfain = polyval(q, n);
    ci(i+1) = alfain / din;
    format rat;
end
end
