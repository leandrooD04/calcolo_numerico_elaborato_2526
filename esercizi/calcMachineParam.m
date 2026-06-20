function [prec, maxNum, minNor, minDen] = calcMachineParam()
% [prec, maxNum, minNor, minDen] = calcMachineParam()
% Funzione che calcola valori base di una aritmetica finita in base binaria
% con arrotondamento
% Output:
%   prec - precisione di macchina
%   maxNum - numero di macchina più grande
%   minNor - numero di macchina più piccolo normalizzato
%   minDen - numero di macchina più piccolo denormalizzato

% 1. Calcolo della precisione
% https://dispense.dmsa.unipd.it/mazzia/prec_e_instabilita.pdf
prec = 1;
while (1.0 + prec / 2) > 1.0
    prec = prec / 2;
end

% 2. Calcolo del numero di macchina più grande
maxNum = 1.0;
while isfinite(maxNum * 2)
    maxNum = maxNum * 2;
end
maxNum = maxNum * (2.0 - prec);

% 3. Calcolo del numero di macchina più piccolo normalizzato
minNor = 1.0;
while (minNor/2.0) + (minNor/2.0) * prec > (minNor/2.0)
    minNor = minNor / 2.0;
end

% 4. Calcolo del numero di macchina più piccolo denormalizzato
minDen = minNor;
while (minDen/2.0) > 0
    minDen = minDen / 2.0;
end

fprintf('Validazione Risultati:\n');
fprintf('precisione: calcolato = %e, eps = %e\n', prec, eps);
fprintf('maxNum: calcolato = %e, realmax = %e\n', maxNum, realmax);
fprintf('minNor: calcolato = %e, realmin = %e\n', minNor, realmin);
fprintf('minDen: calcolato = %e, realmin*eps = %e\n', minDen, realmin*eps);
end