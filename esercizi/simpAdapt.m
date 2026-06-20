function I2 = simpAdapt(f, a, b, tol, fa, fb, fc)
% I2 = simpAdapt(f, a, b, tol, fa, fb, fc)
% Metodo di calcolo dell'integrale di una funzione f per intervallo [a, b]
% tramite la formula adattiva di Simpson
if nargin < 4, error('Numero parametri insufficiente'), end
if nargin < 7, fc = feval(f, (a+b)/2); end
if nargin < 6, fb = feval(f, b); end
if nargin < 5, fa = feval(f, a); end
if a >= b, error('Estremi dell''intervallo non validi'), end
if tol <= 0, error('Tolleranza non valida'), end

h = (b-a)/2;
I1 = (h/3) * (fa+4*fc+fb);

c = (a+b)/2;
d = (a+c)/2;
e = (c+b)/2;
fd = feval(f, d);
fe = feval(f, e);

I2 = (h/6) * (fa + 4*fd + 2*fc + 4*fe + fb);
err = abs(I2 - I1)/15;

if err > tol
    I2 = simpAdapt(f, a, c, tol / 2, fa, fc, fd) + ...
        simpAdapt(f, c, b, tol / 2, fc, fb, fe);
end
end