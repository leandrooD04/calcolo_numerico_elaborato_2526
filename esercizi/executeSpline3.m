function executeSpline3()
% executeSpline3()
% Metodo di esecuzione della function spline3

a = -10;  
b = 10;
xx  = linspace(a, b, 10001);
fxx = runge(xx);

ns = 100:100:2000;
h  = 20 ./ ns;

errNat  = zeros(size(ns));
errKnot = zeros(size(ns));
errComp = zeros(size(ns));

for i = 1:length(ns)
    n  = ns(i);
    xi = linspace(a, b, n+1);
    fi = runge(xi);

    yNat  = spline3(xi, fi, xx, 0); % naturale
    yKnot = spline3(xi, fi, xx, 1);% not-a-knot
    yComp = spline3(xi, fi, xx, 2, drunge(a), drunge(b)); % completa

    errNat(i)  = max(abs(yNat  - fxx));
    errKnot(i) = max(abs(yKnot - fxx));
    errComp(i) = max(abs(yComp - fxx));
end

% Grafico loglog
figure(1); clf;
loglog(h, errNat,  '-o', 'Color', [0.30, 0.75, 1.00], 'LineWidth', 1.6, 'MarkerSize', 5); hold on;
loglog(h, errKnot, '-s', 'Color', [1.00 0.60 0.20], 'LineWidth', 1.6, 'MarkerSize', 5);
loglog(h, errComp, '-^', 'Color', [0.40 1.00 0.40], 'LineWidth', 1.6, 'MarkerSize', 5);

set(gca, 'XDir', 'reverse');
xlabel('h');  ylabel('Errore massimo');
title('Errore spline cubica - funzione di Runge (nodi equispaziati)');
legend('Naturale', 'Not-a-knot', 'Completa', ...
       'pendenza h^2', 'pendenza h^4', 'Location', 'northwest');
grid on;
end

% Funzioni ausiliarie
function y  = runge(x),  y  =  1 ./ (1 + x.^2);        end
function dy = drunge(x), dy = -2*x ./ (1 + x.^2).^2;   end