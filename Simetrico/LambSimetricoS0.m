%% =====================================
% Equa�ao de dispersao de Lamb - Modos Simetricos
%
% Entrada:
% Omega - frequencia angular [rad/s]
% c 	- velocidade de fase [m/s]
% mat	- propriedades do material
%
% Saida:
% F	- funcao caracterestica
%% =====================================
function F = lambSimetricoS0(omega,c,mat)
%% =====================================
% Velocidade do material
cL = mat.cL;
cT = mat.cT;
%% =====================================
% Geometria
h = mat.d/2;
%% =====================================
% Numero de onda
k = omega/c;
%% =====================================
% Numeros de onda na direcao da espessura
p = sqrt(complex((omega/cL)^2 - k^2));
q = sqrt(complex((omega/cT)^2 - k^2));
%% =====================================
% Correcao para valores complexos
% Evita problemas numericos quando o termo interno fica negativo
% if imag(p) ~=0
% 	p = -1i*abs(imag(p));
% end
% if imag(q) ~=0
% 	q = -1i*abs(imag(q));
% end
%% ======================================
% Equacao caracterestica simetrica
F = (q^2-k^2)^2 .* cos(p*h).*sin(q*h) + 4*k^2*p*q .* cos(q*h).*sin(p*h);
end