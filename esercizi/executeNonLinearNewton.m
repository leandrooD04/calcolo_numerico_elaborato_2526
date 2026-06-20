function executeNonLinearNewton()
% executeNonLinearNewton()
% Esecuzione della function di risoluzione di sistemi non lineari
% con Newton
X0 = {
    [0.4; 2; 0.4; 2],
    [0.4; -2; -0.4; 2],
    [1.4; -2; -1.4; 2],
    [1.4; 2; 1.4; 2],
};

for k = 1:4
    x0 = X0{k};
    [x, iter] = nonLinearNewton(@systemFunction, x0, 1e-12, 1000);
    fprintf('Soluzione caso %d individuata in %d iterazioni \n', k, iter);
    disp(x);
end

fprintf('Tentativo con il vettore nullo:\n');

try
    x = nonLinearNewton(@systemFunction, zeros(4,1), 1e-12, 1000);
    fprintf('Soluzione caso %d individuata in %d iterazioni \n', k, iter);
    disp(x);
catch ME
    fprintf('Errore durante la risoluzione: %s\n', ME.message);
end
end

function [Fx,Jx] = systemFunction(x)

Fx = [
    x(1)^2 + x(2)^2 - 5;
    x(3)^2 + x(4)^4 - 5;
    cos(pi*x(1)) + cos(pi*x(3));
    x(2)^2 - x(4)^2
];

Jx = [
    2*x(1), 2*x(2), 0, 0;
    0, 0, 2*x(3), 4*x(4)^3;
    -pi*sin(pi*x(1)), 0, -pi*sin(pi*x(3)), 0;
    0, 2*x(2), 0, -2*x(4)
];

end