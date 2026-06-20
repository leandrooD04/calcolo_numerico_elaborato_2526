function x = LDLSolve(A, b)
% x = LDLSolve(A, b)
% Metodo di risoluzione di un sistema lineare con fattorizzazione LDL^T
% Input:
%   A - Matrice simmetrica definita positiva dei coefficienti
%   b - vettore dei termini noti
% Output:
%   x - vettore soluzione del sistema
if nargin < 2, error('Parametri insufficienti.'), end
[m, n] = size(A);
if m ~= n || m ~= length(b), error('Dimensione dei parametri non validi'), end

% Fattorizzazione 
if A(1, 1) <= 0, error('Matrice non sdp'), end
A(2:n, 1) = A(2:n, 1) / A(1, 1);

for i = 2:n
    v = (A(i,1:i-1).') .* diag(A(1:i-1, 1:i-1));
    A(i, i) = A(i, i) - A(i, 1:i-1) * v;
    if A(i,i) <= 0, error('Matrice non simmetrica definita positiva'), end
    A(i+1:n, i) = (A(i+1:n, i) - A(i+1:n, 1:i-1) * v) / A(i,i);
end

% Risoluzione sistema Ax = b
x = b;
for i = 1:n
    x(i+1:n) = x(i+1:n) - A(i+1:n, i) * x(i);
end
x = x ./ diag(A);

for i = n:-1:2
    x(1:i-1) = x(1:i-1) - A(i, 1:i-1).'*x(i);
end
end
