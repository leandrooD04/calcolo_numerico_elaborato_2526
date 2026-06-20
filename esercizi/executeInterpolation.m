function executeInterpolation()
% executeInterpolation()
% Esegue calcoli di polinomi interpolanti con metodo di Lagrange, Newton e
% Hermite

a    = -10;
b    = 10;
xx   = linspace(a, b, 10001);
fxx  = runge(xx);
f1xx = drunge(xx);

nsLag  = 10:10:140;
nsHerm = 10:10:70;

errLag = zeros(size(nsLag));
errNew = zeros(size(nsLag));
errH   = zeros(size(nsHerm));
errHy1 = zeros(size(nsHerm));

for i = 1:length(nsLag)
    n  = nsLag(i);
    xi = linspace(a, b, n+1);
    fi = runge(xi);
    errLag(i) = max(abs(lagrange(xi, fi, xx) - fxx));
    errNew(i) = max(abs(newton(xi,  fi, xx) - fxx));
end

for j = 1:length(nsHerm)
    n   = nsHerm(j);
    xi  = linspace(a, b, n+1);
    fi  = runge(xi);
    fi1 = drunge(xi);
    [yh, y1h] = hermite(xi, fi, fi1, xx);
    errH(j)   = max(abs(yh  - fxx));
    errHy1(j) = max(abs(y1h - f1xx));
end

figure(1); clf;
semilogy(nsLag,  errLag, '-o',  'Color', [0.30 0.75 1.00], 'LineWidth', 1.8); hold on; % azzurro
semilogy(nsLag,  errNew, '--s', 'Color', [1.00 0.60 0.20], 'LineWidth', 1.8);          % arancio
semilogy(nsHerm, errH,   '-^',  'Color', [0.40 1.00 0.40], 'LineWidth', 1.8);          % verde
semilogy(nsHerm, errHy1, '--d', 'Color', [1.00 1.00 0.30], 'LineWidth', 1.8);          % giallo
xlabel('n');
ylabel('Errore massimo');
title('Errore interpolazione – funzione di Runge (nodi equispaziati)');
legend('Lagrange', 'Newton', 'Hermite su f', 'Hermite su f''', ...
       'Location', 'northwest');
grid on;
xticks(nsLag);
end

function y  = runge(x),  y  =  1 ./ (1 + x.^2);        end
function dy = drunge(x), dy = -2*x ./ (1 + x.^2).^2;   end