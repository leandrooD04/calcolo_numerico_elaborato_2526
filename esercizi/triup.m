function x = triup(U, b)
% x = triup(U, b)
% Metodo di risoluzione di un sistema triangolare superiore
% Input:
%   U - vettore rappresentante la matrice triangolare superiore dei coefficienti
%   b - vettore dei termini noti
% Output:
%   x - vettore soluzione del sistema

if nargin < 2, error('Numero parametri insufficiente'), end
n = (-1 + sqrt(1 + 8 * length(U))) / 2;
if mod(n, 1) ~= 0, error('Lunghezza vettore della matrice dei coefficienti non valida'), end
n = round(n);
if n ~= length(b), error('Le dimensioni dei vettori non sono valide'), end

U = U(:);   % ← forza U a colonna
x = b(:);
for i = n:-1:1
    idxDiag = i*(i+1)/2;
    if U(idxDiag) == 0, error('Matrice singolare'), end
    x(i) = x(i) / U(idxDiag);
    if i > 1
        idxStart = (i-1)*i/2 + 1;
        idxEnd   = idxDiag - 1;
        x(1:i-1) = x(1:i-1) - U(idxStart:idxEnd) * x(i);
    end
end
end