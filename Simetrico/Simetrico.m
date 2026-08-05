clear all
close all
clc
%% ===============================
% Dispersão ondas de Lamb
%
% Bruno Grebin Gasparotto
% Julho/Agosto 2026
%
% Implementação cálculo curvas
%% ===============================
clear
close all
%% ===============================
% Material
E = 210e9; %Pa
nu = 0.30; %Poisson
rho = 7800; %kg/m^3
%% ===============================
% Geometria
d = 1e-3; %espessura placa (m)
h = d/2;
%% ===============================
%% Propriedades Elásticas
mat.lambda 	= E*nu/((1+nu)*(1-2*nu));
mat.mu		= E/(2*(1+nu));
mat.cL 		= sqrt((mat.lambda+2*mat.mu)/rho);
mat.cT		= sqrt(mat.mu/rho);
%% ===============================
% Velocidade extensional
CE = sqrt(E/(rho*(1-nu^2)));
%% ===============================
% Velocidade de Rayleigh
syms x
eq_R = x^6 ...
    - 8*x^4 ...
    + (24 - 16*(mat.cT/mat.cL)^2)*x^2 ...
    - 16*(1-(mat.cT/mat.cL)^2);
sol = vpasolve(eq_R,x,[0 1]);
CR = double(sol)*mat.cT;
fprintf('Velocidade Rayleigh CR = %.2f m/s\n',CR)
fprintf('Velocidade extensional CE = %.2f m/s\n',CE)
fprintf('Propriedades do material \n');
fprintf('Velocidade longitudinal cL = %.2f m/s \n',mat.cL);
fprintf('Velocidade tangencial   cT = %.2f m/s \n',mat.cT);
%% ===============================
% Análise em frequência
% faixa do produto frequência x espessura
% unidade Hz.mm
fd_min =    0.001;
fd_max =   40.000;
Nf     = 8000    ;
fd = linspace(fd_min,fd_max,Nf);
% conversão para Hz
f = fd.*1e6/(d*1e3);
omega = 2*pi*f;
fprintf('\n');
fprintf('Análise de frequência \n');
fprintf('fd mínimo = %.2f MHz.mm \n',fd_min);
fprintf('fd máximo = %.2f MHz.mm \n',fd_max);
fprintf('f = %.2f Hz\n',f(2000));
%% ================================
% Análise velocidade de fase
c_min 	= 2000; % m/s
c_max 	= 7000; % m/s
Nc 	= 5000;
c = linspace(c_min,c_max,Nc);
fprintf('\n');
fprintf('Análise de velocidade de fase \n');
fprintf('c mínimo =  %.2f m/s \n',c_min);
fprintf('c máximo = %.2f m/s \n',c_max);
%% ===============================
% Resolução da Equação de Lamb
F = zeros(Nf,Nc);
for i = 1:Nf
    for j = 1:Nc
        k = omega(i)/c(j);
        p = sqrt((omega(i)/mat.cL)^2 - k^2);
        q = sqrt((omega(i)/mat.cT)^2 - k^2);
        % Correção para valores imaginários
        if imag(p) ~= 0
            p = -1i*abs(imag(p));
        end
        if imag(q) ~= 0
            q = -1i*abs(imag(q));
        end
        F(i,j) = (q^2-k^2)^2*cos(p*h)*sin(q*h) ...
               + 4*k^2*p*q*cos(q*h)*sin(p*h);
    end
end
%% ===============================
% Verificação da função em uma frequência escolhida
fd_verifica = 5.6182;      % MHz.mm
[~,i_verifica] = min(abs(fd-fd_verifica));
Fteste = F(i_verifica,:);   % todos os valores de F para essa frequência
%% Figura Re(F)
figure
plot(c,real(Fteste),'LineWidth',1.5)
grid on
xlabel('Velocidade de fase (m/s)')
ylabel('Re(F)')
title(sprintf('Parte real - fd = %.4f MHz.mm',fd(i_verifica)))
%% Figura Im(F)
figure
plot(c,imag(Fteste),'LineWidth',1.5)
grid on
xlabel('Velocidade de fase (m/s)')
ylabel('Im(F)')
title(sprintf('Parte imaginária - fd = %.4f MHz.mm',fd(i_verifica)))
%% Figura |F|
figure
plot(c,abs(Fteste),'LineWidth',1.5)
grid on
xlim([3000 7000])
title(sprintf('|F| - fd = %.5f MHz.mm',fd(i_verifica)))
xlabel('Velocidade de fase (m/s)')
ylabel('|F|')
%% ===============================
% Busca das raízes para todas as frequências
maxRaizes = 20;               % escolha um valor suficientemente grande
raizes = NaN(Nf,maxRaizes);
for i = 1:Nf
    % Mistura Im(F) abaixo de cT e Re(F) acima de cT
    Fbusca = real(F(i,:));
    idx = c < mat.cT;
    Fbusca(idx) = imag(F(i,idx));
    % Detecta mudanças de sinal
    idx_raiz = find(Fbusca(1:end-1).*Fbusca(2:end) < 0);
    % Velocidade aproximada das raízes
    c_raiz = (c(idx_raiz)+c(idx_raiz+1))/2;
    % Remove raiz em cT
    tol = 0.5;
    c_raiz(abs(c_raiz - mat.cT) < tol) = [];
    % Armazena na matriz
    nRaizes = length(c_raiz);
    if nRaizes > maxRaizes
        warning('Número de raízes maior que maxRaizes.');
        nRaizes = maxRaizes;
    end
    raizes(i,1:nRaizes) = c_raiz;
end
%% ===============================
% Ramo S0
c_S0 = raizes(:,1);   % velocidade de fase do modo S0
%% ===============================
% Curva de dispersão S0
figure
plot(fd,c_S0,'LineWidth',1.5)
hold on
yline(CE,'--','C_E extensional')
yline(CR,'--','C_R Rayleigh')
grid on
xlabel('fd (MHz.mm)')
ylabel('Velocidade de fase (m/s)')
title('Dispersão')
legend
%% ===============================
% Curvas de dispersão S0 até S5
figure
hold on
% Curvas dos modos
for modo = 1:6
    plot(fd,raizes(:,modo),'LineWidth',1.5)
end
% Linhas de corte teóricas
h_corte = gobjects(5,1);
for n = 1:5
    fd_corte = n*mat.cT/1e3;
    h_corte(n) = xline(fd_corte,'--','LineWidth',1.2);
end
for n = 1:5
    xline(n*mat.cT/1e3,'--',sprintf('fd_{S%d}',n))
end
% Referências horizontais
h_CE = yline(CE,'--');
h_CR = yline(CR,'--');
grid on
xlabel('fd (MHz.mm)')
ylabel('Velocidade de fase (m/s)')
title('Dispersão - Modos Simétricos')
legend([...
    plot(nan,nan,'-'),...
    plot(nan,nan,'-'),...
    plot(nan,nan,'-'),...
    plot(nan,nan,'-'),...
    plot(nan,nan,'-'),...
    plot(nan,nan,'-'),...
    h_corte(1),...
    h_corte(2),...
    h_corte(3),...
    h_corte(4),...
    h_corte(5),...
    h_CE,...
    h_CR],...
    {'S0','S1','S2','S3','S4','S5',...
     'Corte S1','Corte S2','Corte S3','Corte S4','Corte S5',...
     'C_E','C_R'})