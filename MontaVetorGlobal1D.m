function F = MontaVetorGlobal1D(x,u0,ul,ff)
n = length(x)-1; % n é o número de subintervalos - elementos
F = zeros(n+1,1); % aloca vetor
for i = 1:n % loop sobre os elementos
h = x(i+1) - x(i); %tamanho do elemento
F(i) = F(i) + ff(x(i))*h/2;
F(i+1) = F(i+1) + ff(x(i+1))*h/2;
end
F(1) = u0;
h = x(2) - x(1);
F(2) = F(2) - (-1/h)*u0;
F(n+1) = ul;
h = x(n+1) - x(n);
F(n) = F(n) - (-1/h)*ul;
return
