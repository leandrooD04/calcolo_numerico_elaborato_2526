function [x, nor] = leastSquareSolve(A, b, omega)
% x = leastSquareSolve(A, b, omega)
% Metodo di risoluzione di un sistema lineare sovradeterminato con
% soluzione ai minimi quadrati pesati
% Input:
%   A - matrice sovradeterminata dei coefficienti
%   b - vettore dei termini noti
%   omega - vettore dei pesi positivi
% Output:
%   x - vettore soluzione del sistema
%   nor - norma del vettore residuo
if nargin < 3, error('Numero parametri insufficiente'), end
[m, n] = size(A);
if m <= n, error('Dimensione matrice non valida'), end
if m ~= length(b) || m ~= length(omega), error('Dimensione del vettore dei termini noti non valida'), end
if any(omega <= 0), error('I pesi non sono tutti positivi'), end

W = sqrt(omega(:));
A = A .* W;
b = b .* W;

for i = 1:n
    alpha = norm(A(i:m, i));
    if alpha == 0, error('Matrice non ha il rango massimo'), end
    if A(i,i) >= 0, alpha = -alpha; end
    v1 = A(i,i) - alpha;
    A(i,i) = alpha;
    A(i+1:m, i) = A(i+1:m, i)/v1;
    beta = -v1/alpha;
    
    v = [1;A(i+1:m, i)];
    % applica H = I - beta*v*v' alle colonne successive
    A(i:m, i+1:n) = A(i:m, i+1:n) - (beta * v * (v' * A(i:m, i+1:n)));

    % applica lo stesso riflettore a b
    b(i:m) = b(i:m) - (beta * v * (v' * b(i:m)));
end

% Sostituzione all'indietro:
x = zeros(n, 1);
for i = n:-1:1
    x(i) = (b(i) - A(i, i+1:n) * x(i+1:n)) / A(i,i);
end
nor = norm(b(n+1:m));
end