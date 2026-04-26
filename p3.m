clear all
clc

% Primeiro define n, depois cria o vetor orientation (corrigido)
n = input('Insira o número de camadas (inteiro): '); % Número de camadas que o usuário deve definir
t = input('Insira a espessura das camadas: '); % Tamanho de espessura t que todas as camadas devem ter
orientation = zeros(1, n); % Corrigido: agora n já está definido

% Preenche as orientações de cada camada (corrigido para for e até n)
for i = 1:n
    mensagem = sprintf('Insira a orientação (em graus) da camada %d: ', i); % Graus
    orientation(i) = input(mensagem);
end

E1  = input('Insira o  modelo de elasticidade E1 ');  %Pa
E2  = input('Insira o  modelo de elasticidade E2 ');  %Pa

v12 = input('Insira o  coeficiente de Poisson v12 '); %Pa
v13 = input('Insira o  coeficiente de Poisson v13 '); %Pa
v23 = input('Insira o  coeficiente de Poisson v23 '); %Pa

G12 = input('Insira o coeficiente de cisalhamento G12 ');  %Pa
G13 = input('Insira o coeficiente de cisalhamento G13 ');  %Pa
G23 = input('Insira o coeficiente de cisalhamento G23 ');  %Pa

alpha1  = input('Insira o coeficiente de expansão térmica 1 ');  %Pa
alpha2  = input('Insira o coeficiente de expansão térmica 2 ');  %Pa

Xet = input('Insira a resistência da lâmina Xet ');  %Pa
Xec = input('Insira a resistência da lâmina Xec ');  %Pa
Yet = input('Insira a resistência da lâmina Yet ');  %Pa
Yec = input('Insira a resistência da lâmina Yec ');  %Pa
Se  = input('Insira a resistência da lâmina Se ');   %Pa

Q = zeros(3, 3, n); % Matriz de rigidez de cada camada
alphaXY = zeros(3, 1, n); % Corrigido: agora armazena 3 componentes para cada camada. Direção X, Y e XY

v21 = v12 * E2 / E1; % Corrigido: cálculo de v21

Q11 = E1/(1 - v12*v21);
Q12 = v12*E2/(1 - v12*v21);
Q22 = E2/(1 - v12*v21);
Q66 = G12;

for i = 1:n
    theta = deg2rad(orientation(i)); % Corrigido: variável theta para não sobrescrever t
    c = cos(theta); % Cosseno da orientação
    s = sin(theta); % Seno da orientação

    % Matriz de transformação para as tensões
    T = [c^2 s^2 2*c*s; s^2 c^2 -2*c*s; -c*s c*s c^2-s^2]; 

    % Cálculo dos elementos da matriz Q para cada camada (corrigido para Q(:,:,i))
    Q(1,1,i) = Q11*c^4 + 2*(Q12 + 2*Q66)*c^2*s^2 + Q22*s^4; % Elemento Q11
    Q(1,2,i) = (Q11 + Q22 - 4*Q66)*c^2*s^2 + Q12*(c^4 + s^4); % Elemento Q12
    Q(1,3,i) = (Q11 - Q12 - 2*Q66)*c^3*s - (Q22 - Q12 - 2*Q66)*c*s^3; % Elemento Q13
    Q(2,2,i) = Q11*s^4 + 2*(Q12 + 2*Q66)*c^2*s^2 + Q22*c^4; % Elemento Q22
    Q(2,3,i) = (Q11 - Q12 - 2*Q66)*c*s^3 - (Q22 - Q12 - 2*Q66)*c^3*s; % Elemento Q23
    Q(3,3,i) = (Q11 + Q22 - 2*Q12 - 2*Q66)*c^2*s^2 + Q66*(c^4 + s^4); % Elemento Q33

    % Demonstração da simetria da matriz
    Q(2,1,i) = Q(1,2,i); 
    Q(3,1,i) = Q(1,3,i); 
    Q(3,2,i) = Q(2,3,i);

    % Cálculo do vetor alphaXY (corrigido para alphaXY(3,1,n))
    alphaXY(1,1,i) = alpha1*c^2 + alpha2*s^2; % Elemento X de alpha
    alphaXY(2,1,i) = alpha1*s^2 + alpha2*c^2; % Elemento Y de alpha
    alphaXY(3,1,i) = (alpha1 - alpha2)*c*s; % Elemento XY de alpha
end

Nx  = input('Insira o carregamento mecânico em x:');   %Pa
Ny  = input('Insira o carregamento mecânico em y:');   %Pa
Nxy = input('Insira o carregamento mecânico em xy:');   %Pa

Mx  = input('Insira o carregamento mecânico em x:');   %Pa
My  = input('Insira o carregamento mecânico em y:');   %Pa
Mxy = input('Insira o carregamento mecânico em xy:');   %Pa

DT = input('Insira a variação de temperatura DT (em graus)');   %graus Celsius

altura = (n*t); % Altura total do laminado

z = zeros(1, n+1);
z(1) = -altura/2;
for i = 1:n
    z(i+1) = z(i)+ t;
end

% Inicialização das matrizes 
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
As = zeros(2,2);

% Removido input duplicado de G13 e G23

for i = 1:n
    % Matrizes A, B e D (corrigido para Q(:,:,i))
    A = A + Q(:,:,i) * (z(i+1) - z(i));
    B = B + 0.5 * Q(:,:,i) * (z(i+1)^2 - z(i)^2);
    D = D + (1/3) * Q(:,:,i) * (z(i+1)^3 - z(i)^3);

    % Matriz As (cisalhamento transversal)
    theta = deg2rad(orientation(i)); % Corrigido: variável theta para não sobrescrever t
    c = cos(theta); % Cosseno da orientação
    s = sin(theta); % Seno da orientação

    % Matriz de rigidez de cisalhamento local
    Qs = [G13, 0; 0 , G23];

    % Matriz de rotação para o sistema local
    Ts = [c, s; -s, c];
    Qbs = Ts.' * Qs * Ts;

    As = As + Qbs * (z(i+1) - z(i));
end


% =========================
% Etapa 8: Cálculo das deformações na superfície média (e0) e curvaturas (kapa)
% =========================

% Vetores de carregamento
N = [Nx; Ny; Nxy]; % Forças normais e de cisalhamento
M = [Mx; My; Mxy]; % Momentos

% Monta a matriz ABD
ABD = [A B; B D];

% Vetor de carregamento total
NM = [N; M];

% Resolve o sistema para encontrar e0 e kapa
sol = ABD \ NM;
e0 = sol(1:3); % Deformações na superfície média (e0x, e0y, Y0xy)
kapa = sol(4:6); % Curvaturas (kapa_x, kapa_y, kapa_xy)

% Exibe resultados
disp('Deformações na superfície média (e0):');
disp(e0);
disp('Curvaturas (kapa):');
disp(kapa);

% =========================
% Etapas 9 e 10: Definição de npl e cálculo das posições z dos pontos analisados
% =========================

% Solicita ao usuário o número de pontos por camada (npl >= 3)
npl = input('Informe a quantidade de pontos por camada para análise (npl >= 3): ');
while npl < 3
    disp('O valor de npl deve ser maior ou igual a 3.');
    npl = input('Informe a quantidade de pontos por camada para análise (npl >= 3): ');
end

% Calcula as posições z dos pontos analisados ao longo do laminado
totalPoints = n * npl;
zpos = zeros(1, totalPoints);

idx = 1;
for camada = 1:n
    z_base = z(camada);
    dz = t / (npl - 1);
    for p = 0:(npl-1)
        zpos(idx) = z_base + p * dz;
        idx = idx + 1;
    end
end
