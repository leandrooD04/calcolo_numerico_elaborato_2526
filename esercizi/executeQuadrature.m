function executeQuadrature()
% executeQuadrature()
% Approssima  int_0^1 100*sin(100x)*exp(-100x) dx  (~ 1/2)  con la formula
% composita di Newton-Cotes e adattiva di Simpson, con tabulazione dei
% dati.

f = @(x) 100*sin(100*x).*exp(-100*x);
a = 0;  
b = 1;  
Iexact = 0.5;

nfeval = 0;
g = @count;

% ---- formula composita di Newton-Cotes ----
fprintf('Formula composita di Newton-Cotes \n');
fprintf('%6s %12s %18s %14s %14s\n', 'n', 'val.funz.', 'I (approx)', 'errore', 'stima err.');
for n = [40 80 160 320]
    nfeval = 0;
    [Ikn, err] = ncotescomp(g, a, b, 2, n);
    fprintf('%6d %12d %18.12f %14.3e %14.3e\n', n, nfeval, Ikn, abs(Ikn - Iexact), abs(err));
end

% ---- formula adattiva di Simpson ----
fprintf('\nFormula adattiva di Simpson\n');
fprintf('%8s %12s %18s %14s\n', 'tol', 'val.funz.', 'I (approx)', 'errore');
for tol = [1e-3 1e-4 1e-5 1e-6]
    nfeval = 0;
    I = simpAdapt(g, a, b, tol);
    fprintf('%8.0e %12d %18.12f %14.3e\n', tol, nfeval, I, abs(I - Iexact));
end

    % funzione annidata: valuta f e aggiorna il contatore di executeQuadrature
    function y = count(x)
        nfeval = nfeval + numel(x);
        y = f(x);
    end
end