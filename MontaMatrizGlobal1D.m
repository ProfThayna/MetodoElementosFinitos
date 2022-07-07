function A = MontaMatrizGlobal1D(x)
n = length(x)-1; % n é o número de subintervalos - elementos
A = zeros(n+1,n+1); % aloca matriz de rigidez
for i = 1:n % loop sobre os elementos
h = x(i+1) - x(i); %tamanho do elemento
A(i,i) = A(i,i) + 1/h;
A(i,i+1) = A(i,i+1) - 1/h;
A(i+1,i) = A(i+1,i) - 1/h;
A(i+1,i+1) = A(i+1,i+1) + 1/h;
end
A(1,1) = 1;
A(n+1,n+1) = 1;
A(1,2) = 0;
A(2,1) = 0;
A(n+1,n) = 0;
A(n,n+1) = 0;
return
