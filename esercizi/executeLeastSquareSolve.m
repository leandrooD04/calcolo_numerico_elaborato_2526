function executeLeastSquareSolve()
% executeLeastSquareSolve()
% Esecuzione della function di risoluzione di sistemi lineari
% sovradeterminati pesati
A = [1 2 3; 4 5 6; 7 8 9; 10 11 0];
b = [1; 2; 3; 4];
i = (1:4)';

omegas = {
    ones(4,1);
    ((5-i).^2)/9;
    (i.^2)/9;
    sqrt(2)*ones(4,1);
};

for k = 1:4
    [x, nor] = leastSquareSolve(A, b, omegas{k});
    fprintf('Caso %d:\n', k);
    fprintf('   pesi    = [%s]\n', sprintf('%.4f ', omegas{k}));
    fprintf('   x       = [%s]\n', sprintf('%.6f ', x));
    fprintf('   ||r||_w = %.6e\n\n', nor);
end
end