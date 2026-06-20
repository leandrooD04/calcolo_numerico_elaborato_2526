function x = LUPivotingSolve(A, b)
% x = LUPivotingSolve(A, b)
% Metodo di risoluzione di sistema lineare tramite fattorizzazione LU con
% pivoting
% Input:
%   A - matrice dei coefficienti
%   b - vettore dei termini noti
% Output:
%   x - vettore soluzione del sistema

if nargin < 2, error('Numero di parametri insufficiente'), end
[m, n] = size(A);
if m ~= n || m ~= length(b), error('Dimensioni dei parametri non validi'), end

P = (1:n);
% Fattorizzazione LU con Pivoting parziale
for i = 1:n-1
    [mi, ki] = max(abs(A(i:n, i)));
    if mi == 0, error('Matrice singolare'), end
    ki = ki + i - 1;
    if ki > i
        A([i ki], :) = A([ki i], :);
        P([i ki]) = P([ki i]);
    end
    A(i+1:n, i) = A(i+1:n, i) / A(i, i);
    A(i+1:n, i+1:n) = A(i+1:n, i+1:n) - A(i+1:n, i) * A(i, i+1:n);
end

if A(n, n) == 0, error('Matrice singolare'), end
x = b(P);

% Risoluzione del sistema per L
for i = 2:n
    x(i:n) = x(i:n) - A(i:n, i-1) * x(i-1);
end

% Risoluzione del sistema per U
for i = n:-1:1
    x(i) = x(i) / A(i, i);
    x(1:i-1) = x(1:i-1) - A(1:i-1, i) * x(i);
end
end