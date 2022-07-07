clear all
clc
close all
L = 1.0; % comprimento do domínio
u0=0; % valor de contorno em x=0
ul=0.0; % valor de contorno em x=l
n = 100; % numero de elementos (subintervalos)
h = L/n; % comprimento dos elementos (iguais)
x = 0:h:L; % malha - vetor de pontos nodais
A = MontaMatrizGlobal1D(x); % monta matriz global
F = MontaVetorGlobal1D(x,u0,ul,@ff); % monta vetor global
ue = sin(pi.*x);
ue = ue';
u = A\F; % resolve o sistema linear
figure,
plot(x,u,'r.',x, ue,'k');
title('Solução MEF');
[u ue];
figure,
plot(x, abs(u-ue))
title('Erro');



