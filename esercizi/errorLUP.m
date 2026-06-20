function [x, err] = errorLUP()
% [x, err] = errorLUP()
% Metodo di verifica dell'errore della fattorizazione LU con pivoting.
% Genera 15 matrici A e vettori b e risolve i relativi sistemi lineari.
% Output:
%   x - vettore soluzione
%   err: errore commesso

n = 10;
y = ones(n, 1);

for k = 1:15
    [A,b] = linsis(n,k);

    x = LUPivotingSolve(A, b);
    err = norm(x-y)/norm(y);
    c = cond(A);
    fprintf('Risultati sistema lineare numero %d: errore = %e, condizionamento = %e\n', k, err, c);
end
end