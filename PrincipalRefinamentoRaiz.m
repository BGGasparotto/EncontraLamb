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
%% ===============================
%% Propriedades Elásticas
mat = MaterialS0 (E,nu,rho,d);
fprintf('Propriedades do material \n');
fprintf('Velocidade longitudinal cL = %.2f m/s \n',mat.cL);
fprintf('Velocidade tangencial   cT = %.2f m/s \n',mat.cT);
%% ===============================
% Análise em frequência
% faixa do produto frequência x espessura
% unidade Hz.mm
fd_min =    0.001;
fd_max =   10.000;
Nf     = 1000    ;
fd = linspace(fd_min,fd_max,Nf);
% conversão para Hz
f = fd.*1e6/(d*1e3);
omega = 2*pi*f;
fprintf('\n');
fprintf('Análise de frequência \n');
fprintf('fd mínimo = %.2f MHz.mm \n',fd_min);
fprintf('fd máximo = %.2f MHz.mm \n',fd_max);
fprintf('f = %.2f Hz\n',f(100));
%% ================================
% Análise velocidade de fase
c_min 	= 2000; % m/s
c_max 	= 7000; % m/s
Nc 	= 3000;
c = linspace(c_min,c_max,Nc);
fprintf('\n');
fprintf('Análise de velocidade de fase \n');
fprintf('c mínimo =  %.2f m/s \n',c_min);
fprintf('c máximo = %.2f m/s \n',c_max);
%% ===============================
%Resolução da Equação de Lamb
zeros_S = cell(Nf,1);
for i = 1:Nf
    Fteste=zeros(size(c));
    for j=1:Nc
        Fteste(j)=LambSimetricoS0(omega(i),c(j),mat);
    end
    % Remove região abaixo de cT
    idx_ct = c > mat.cT*0.7;
    c_busca = c(idx_ct);
    F_busca = real(Fteste(idx_ct));
if abs(fd(i)-9.70979) < 0.005
    % Parte real
    figure
    plot(c_busca,real(Fteste(idx_ct)),'LineWidth',1.5)
    hold on
    yline(0,'k--')
    grid on
    xlim([3000 7000])
    title(sprintf('Re(F) - fd = %.4f',fd(i)))
    xlabel('Velocidade (m/s)')
    ylabel('Re(F)')
    % Módulo
    figure
    plot(c_busca,abs(Fteste(idx_ct)),'LineWidth',1.5)
    grid on
    xlim([3000 7000])
    title(sprintf('|F| - fd = %.4f',fd(i)))
    xlabel('Velocidade (m/s)')
    ylabel('|F|')
end
    % Procura mudanças de sinal
    idx1 = find(F_busca(1:end-1).*F_busca(2:end)<0);
    raizes1 = [];
for k = 1:length(idx1)
    c1 = c_busca(idx1(k));
    c2 = c_busca(idx1(k)+1);
    raiz = fzero(@(x) real(LambSimetricoS0(omega(i),x,mat)), [c1 c2]);
    raizes1(end+1) = raiz;
end
    %% Junta as raízes
    zeros_S{i}=sort([raizes1]);
end
%% ================================
% Tabela de raízes
Nzeros_max = max(cellfun(@length,zeros_S));
Tabela_raizes = NaN(Nf,Nzeros_max);
for i = 1:Nf
    nr = length(zeros_S{i});
    if nr > 0
        Tabela_raizes(i,1:nr)=zeros_S{i};
    end
end
Tabela_raizes = array2table(Tabela_raizes);
Tabela_raizes.fd = fd';
Tabela_raizes = movevars(Tabela_raizes,'fd','Before',1);
disp(Tabela_raizes)
%% ================================
% Zeros encontrados
    fprintf('Zeros em fd = %.3f MHz.mm\n',fd(end))
for n = 1:length(zeros_S{end})
    fprintf('Zero %d: %.2f m/s\n',...
        n,zeros_S{end}(n))
end
%% ===============================
% Rastreamento direto do modo S0
cS0 = NaN(Nf,1);
% Velocidade extensional
cE = sqrt(E/(rho*(1-nu^2)));
%------------------------------------------
% Primeiro ponto
F = zeros(size(c));
for j = 1:Nc
    F(j) = real(LambSimetricoS0(omega(1),c(j),mat));
end
idx = find(F(1:end-1).*F(2:end)<0);
raizes = [];
for k = 1:length(idx)
    c1 = c(idx(k));
    c2 = c(idx(k)+1);
    raiz = fzero(@(x) real(LambSimetricoS0(omega(1),x,mat)),[c1 c2]);
    raizes(end+1)=raiz;
end
[~,i0]=min(abs(raizes-cE));
cS0(1)=raizes(i0);
fprintf('Primeira raiz S0 = %.2f m/s\n',cS0(1));
%------------------------------------------
% Continuação
janela = 400;      % m/s
for i = 2:Nf
    cmin = max(c_min,cS0(i-1)-janela);
    cmax = min(c_max,cS0(i-1)+janela);
    idxJanela = c>=cmin & c<=cmax;
    cLocal = c(idxJanela);
    F = zeros(size(cLocal));
    for j = 1:length(cLocal)
        F(j)=real(LambSimetricoS0(omega(i),cLocal(j),mat));
    end
    idx = find(F(1:end-1).*F(2:end)<0);
    if isempty(idx)
        continue
    end
    candidatos = zeros(length(idx),1);
    for k = 1:length(idx)
        c1 = cLocal(idx(k));
        c2 = cLocal(idx(k)+1);
        candidatos(k)=fzero(@(x) real(LambSimetricoS0(omega(i),x,mat)),...
                            [c1 c2]);
    end
    [~,k]=min(abs(candidatos-cS0(i-1)));
    cS0(i)=candidatos(k);
end

%% ===============================
% Gráfico

figure
plot(fd,cS0,'LineWidth',2)
grid on
xlabel('fd (MHz.mm)')
ylabel('Velocidade de fase (m/s)')
title('Modo S0')