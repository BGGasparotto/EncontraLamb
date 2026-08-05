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
Nf     = 2000    ;
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
Nc 	= 2000;
c = linspace(c_min,c_max,Nc);
fprintf('\n');
fprintf('Análise de velocidade de fase \n');
fprintf('c mínimo =  %.2f m/s \n',c_min);
fprintf('c máximo = %.2f m/s \n',c_max);
%% ===============================
% Resolução da Equação de Lamb
% Continuação do ramo S0
zeros_S = zeros(Nf,1);
%% ===============================
% Verificação da função em uma frequência escolhida
fd_verifica = 8; % MHz.mm
[~,i_verifica] = min(abs(fd-fd_verifica));
Fteste=zeros(size(c));
for j=1:Nc
    Fteste(j)=(LambSimetricoS0(omega(i_verifica),c(j),mat));
end
% Remove região abaixo de cT
idx_ct = c > mat.cT*0.7;
c_busca = c(idx_ct);
F_busca = Fteste(idx_ct);
%% Figura Re(F)
figure
plot(c_busca,real(Fteste(idx_ct)),'LineWidth',1.5)
hold on
yline(0,'k--')
grid on
xlim([3000 7000])
title(sprintf('Re(F) - fd = %.5f MHz.mm',fd(i_verifica)))
xlabel('Velocidade de fase (m/s)')
ylabel('Re(F)')
%% Figura Im(F)
figure
plot(c_busca,imag(Fteste(idx_ct)),'LineWidth',1.5)
grid on
xlim([3000 7000])
title(sprintf('Im(F) - fd = %.5f MHz.mm',fd(i_verifica)))
xlabel('Velocidade de fase (m/s)')
ylabel('Im(F)')
%% Figura |F|
figure
plot(c_busca,abs(Fteste(idx_ct)),'LineWidth',1.5)
grid on
xlim([3000 7000])
title(sprintf('|F| - fd = %.5f MHz.mm',fd(i_verifica)))
xlabel('Velocidade de fase (m/s)')
ylabel('|F|')
%% ===============================
% Inicialização do ramo S0
% Primeira frequência
i = 1;
Fteste=zeros(size(c));
for j=1:Nc
    Fteste(j)=real(LambSimetricoS0(omega(i),c(j),mat));
end
% Remove região abaixo de cT
idx_ct = c > mat.cT*0.7;
c_busca = c(idx_ct);
F_busca = Fteste(idx_ct);
% Procura mudança de sinal
idx=find(F_busca(1:end-1).*F_busca(2:end)<0);
if isempty(idx)
    error('Nenhuma raiz inicial encontrada para o modo S0');
end
% Primeira raiz do ramo S0
cAnt=fzero(@(x) real(LambSimetricoS0(omega(i),x,mat)),...
           [c_busca(idx(1)) c_busca(idx(1)+1)]);
zeros_S(i)=cAnt;
%% ===============================
% Continuação do ramo
for i = 2:Nf
    % Previsão pela solução anterior
    cPred = cAnt;
    % Janela normal de busca
    dc = 50;
    cmin = cPred - dc;
    cmax = cPred + dc;
    % Limites físicos
    cmin = max(cmin,c_min);
    cmax = min(cmax,c_max);
   %% Avalia função dentro da janela
cLocal = linspace(cmin,cmax,1000);
FLocal = zeros(size(cLocal));
ImLocal = zeros(size(cLocal));
for j=1:length(cLocal)
    Fvalor = LambSimetricoS0(omega(i),cLocal(j),mat);
    FLocal(j)=real(Fvalor);
    ImLocal(j)=imag(Fvalor);
end
 % Procura mudança de sinal
    idx=find(FLocal(1:end-1).*FLocal(2:end)<0);
%% ===============================
% Existe raiz real na janela
if ~isempty(idx)
    cNovo=fzero(...
        @(x) real(LambSimetricoS0(omega(i),x,mat)),...
        [cLocal(idx(1)) cLocal(idx(1)+1)]);
%% ===============================
% Não existe raiz em Re(F)
else
%% Região próxima de cT
if abs(cPred-mat.cT)<=dc
    tol_im = 1e-8;
    idx_zero = find(abs(ImLocal)<=tol_im,1,'first');
    if ~isempty(idx_zero)
        cNovo = cLocal(idx_zero);
    else
        % tenta Re(F)
        if ~isempty(idx)
            cNovo=fzero(...
            @(x) real(LambSimetricoS0(omega(i),x,mat)),...
            [cLocal(idx(1)) cLocal(idx(1)+1)]);
        else
            cNovo=cPred;
        end
    end
else
    % região normal
    if ~isempty(idx)
        cNovo=fzero(...
        @(x) real(LambSimetricoS0(omega(i),x,mat)),...
        [cLocal(idx(1)) cLocal(idx(1)+1)]);
    else
        cNovo=cPred;
    end
end
end
    % Atualiza solução
    zeros_S(i)=cNovo;
    cAnt=cNovo;
    end
%% ================================
% Tabela de raízes S0
Tabela_raizes = table(fd', zeros_S,...
    'VariableNames', {'fd','c_S0'});
Tabela_raizes = movevars(Tabela_raizes,'fd','Before',1);
disp(Tabela_raizes)
%% ================================
%% ================================
% Zero final encontrado
fprintf('Zero S0 em fd = %.3f MHz.mm\n',fd(end))
fprintf('c_S0 = %.2f m/s\n',zeros_S(end))
%% ===============================
figure
plot(fd,zeros_S,'LineWidth',1.5)
grid on
xlabel('fd (MHz.mm)')
ylabel('Velocidade de fase (m/s)')
title('Modo S0')