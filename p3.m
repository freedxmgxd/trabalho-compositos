clear all
clc

% 1 – Definir o número de n camadas do laminado, onde n deve definido pelo programador;

% n = input('Insira o número de camadas (inteiro): '); % Número de camadas que o usuário deve definir
n = 3

% 2 – Todas as camadas devem ter a mesma espessura t que será definida pelo programador;

% t = input('Insira a espessura das camadas: '); % Tamanho de espessura t que todas as camadas devem ter
t = 5.0000e-04

% 3 – Deve ser criado um vetor com as orientações de cada camada seguindo o padrão:
% orientation(1) = ##;
% orientation(2) = ##;
% ...
% orientation(n) = ##;
% onde os valores das orientações devem ser fornecidos em graus pelo programador

orientation = zeros(1, n);

% Preenche as orientações de cada camada
% for i = 1:n
%     mensagem = sprintf('Insira a orientação (em graus) da camada %d: ', i); % Graus
%     orientation(i) = input(mensagem);
% end

orientation = [0, 30, -45]; % Exemplo de orientações para 3 camadas

% 4 – Receber as propriedades da lâmina nas direções 1 e 2 (módulos de elasticidade - E1, E2, coeficientes de Poisson - v12, v13, v23, módulos de cisalhamento - G12, G13, G23, coeficientes de expansão térmica alpha1 e alpha2, resistências da lâmina - Xt,Xc,Yt,Yc,S6 ou Xet,Xec,Yet,Yec,Se);

% E1 = input('Insira o  modelo de elasticidade E1 '); %Pa
% E2 = input('Insira o  modelo de elasticidade E2 '); %Pa

% v12 = input('Insira o  coeficiente de Poisson v12 '); %Pa
% v13 = input('Insira o  coeficiente de Poisson v13 '); %Pa
% v23 = input('Insira o  coeficiente de Poisson v23 '); %Pa

% G12 = input('Insira o coeficiente de cisalhamento G12 '); %Pa
% G13 = input('Insira o coeficiente de cisalhamento G13 '); %Pa
% G23 = input('Insira o coeficiente de cisalhamento G23 '); %Pa

% alpha1 = input('Insira o coeficiente de expansão térmica 1 '); %Pa
% alpha2 = input('Insira o coeficiente de expansão térmica 2 '); %Pa

% Xet = input('Insira a resistência da lâmina Xet '); %Pa
% Xec = input('Insira a resistência da lâmina Xec '); %Pa
% Yet = input('Insira a resistência da lâmina Yet '); %Pa
% Yec = input('Insira a resistência da lâmina Yec '); %Pa
% Se = input('Insira a resistência da lâmina Se '); %Pa

E1 = 1.8100e+11
E2 = 1.0300e+10
% E3 = 1.0300e+10
v12 = 0.2800
v13 = 0.2800
v23 = 0.2800
G12 = 7.1700e+09
G13 = 7.1700e+09
G23 = 7.1700e+09
alpha1 = 0
alpha2 = 0
Xet = 1.0350e+09
Xec = 1.0350e+09
Yet = 2.7600e+07
Yec = 1.3800e+08
Se = 4.1400e+07

% 5 – Para cada camada do laminado calcular as matrizes [𝑄] e o vetor {𝛼}𝑥𝑦.
% Para a matriz [𝑄] usar uma matriz Q=zeros(3,3,n) e para o vetor {𝛼}𝑥𝑦 usar uma matriz alphaXY=zeros(3,1,n), onde:
% Q=(3,3,1) e alphaXY=(3,1,1) são relativos à primeira camada;
% Q=(3,3,2) e alphaXY=(3,1,2) são relativos à segunda camada;
% ...
% Q=(3,3,n) e alphaXY=(3,1,n) são relativos à camada n;

Q = zeros(3, 3, n);
L = zeros(3, 3, n);

alphaXY = zeros(3, 1, n);

v21 = v12 * E2 / E1;

Q0 = [E1 / (1 - v12 * v21), v12 * E2 / (1 - v12 * v21), 0;
    v12 * E2 / (1 - v12 * v21), E2 / (1 - v12 * v21), 0;
    0, 0, G12]

for i = 1:n
    theta = orientation(i);

    c = cosd(theta); % Cosseno da orientação
    s = sind(theta); % Seno da orientação

    L(:, :, i) = [c^2, s^2, c * s;
            s^2, c^2, -c * s;
            -2 * c * s, 2 * c * s, c^2 - s^2];

    Q(:, :, i) = L(:, :, i)' * Q0 * L(:, :, i)

    % Cálculo do vetor alphaXY
    alphaXY(1, 1, i) = alpha1 * c^2 + alpha2 * s^2; % Elemento X de alpha
    alphaXY(2, 1, i) = alpha1 * s^2 + alpha2 * c^2; % Elemento Y de alpha
    alphaXY(3, 1, i) = (alpha1 - alpha2) * c * s; % Elemento XY de alpha
end

% 6 – Receber o estado os carregamentos mecânicos Nmec = (Nx, Ny e Nxy), Mmec = (Mx, My e Mxy) e a variação de temperatura DT em graus Celsius;

% Nx = input('Insira o carregamento mecânico em x:'); %Pa
% Ny = input('Insira o carregamento mecânico em y:'); %Pa
% Nxy = input('Insira o carregamento mecânico em xy:'); %Pa

% Mx = input('Insira o carregamento mecânico em x:'); %Pa
% My = input('Insira o carregamento mecânico em y:'); %Pa
% Mxy = input('Insira o carregamento mecânico em xy:'); %Pa

% % Vetores de carregamento
% Nmec = [Nx; Ny; Nxy]; % Forças normais e de cisalhamento
% Mmec = [Mx; My; Mxy]; % Momentos

Nmec = [150000; 0; 0]

Mmec = [0; 0; 0]

% DT = input('Insira a variação de temperatura DT (em graus)'); %graus Celsius
DT = 0;

% 7 – Calcular as matrizes [A], [B], [D] e [As];
% Utilizando a Teoria Clássica da Laminação:

altura = (n * t); % Altura total do laminado

z = zeros(1, n + 1);
z(1) = altura / 2;

for i = 1:n
    z(i + 1) = z(i) - t;
end

z

% Inicialização das matrizes
A = zeros(3, 3);
B = zeros(3, 3);
D = zeros(3, 3);
As = zeros(2, 2);

for i = 1:n
    % Matrizes A, B e D
    A = A + Q(:, :, i) * (z(i) - z(i + 1));
    B = B + (1/2) * Q(:, :, i) * (z(i)^2 - z(i + 1)^2);
    D = D + (1/3) * Q(:, :, i) * (z(i)^3 - z(i + 1)^3);

    % Matriz As (cisalhamento transversal)
    theta = orientation(i);
    c = cosd(theta); % Cosseno da orientação
    s = sind(theta); % Seno da orientação

    % Matriz de rigidez de cisalhamento local
    Qs = [G13, 0; 0, G23];

    % Matriz de rotação para o sistema local
    Ts = [c, s; -s, c];
    Qbs = Ts.' * Qs * Ts;

    As = As + Qbs * (z(i + 1) - z(i));
end

% 8 – Calcular as deformações na superfície média do laminado:
% e0 = (e0x, e0y, Y0xy);
% kapa = (kapa_x, kapa _y, kapa _xy);

% Monta a matriz ABD
ABD = [A B; B D];

% Vetor de carregamento total
NM = [Nmec; Mmec];

% Resolve o sistema para encontrar e0 e kapa
sol = ABD \ NM;
e0 = sol(1:3); % Deformações na superfície média (e0x, e0y, Y0xy)
kapa = sol(4:6); % Curvaturas (kapa_x, kapa_y, kapa_xy)

% Exibe resultados
fprintf('Deformações na superfície média (e0):\n');
fprintf('  %e\n', e0);
fprintf('Curvaturas (kapa):\n');
fprintf('  %e\n', kapa);

% 9 – Criar a variável npl que define a quantidade de pontos dentro de uma camada onde as deformações, tensões e critérios de falhas serão analisados. O código deve funcionar para qualquer valor de npl. No mínimo usar npl igual a 3.

% Solicita ao usuário o número de pontos por camada (npl >= 3)
% npl = input('Informe a quantidade de pontos por camada para análise (npl >= 3): ');
npl = 3;

while npl < 3
    fprintf('O valor de npl deve ser maior ou igual a 3.');
    npl = input('Informe a quantidade de pontos por camada para análise (npl >= 3): ');
end

% 10 – Calcular o valor da coordenada z para cada ponto analisado ao longo do laminado e armazenar esses valores na variável zpos;

% Calcula as posições z dos pontos analisados ao longo do laminado
totalPoints = n * npl;
zpos = zeros(1, totalPoints);

idx = 1;

for camada = 1:n
    z_base = z(camada);
    dz = t / (npl - 1);

    for p = 0:(npl - 1)
        zpos(idx) = z_base - p * dz;
        idx = idx + 1;
    end

end

zpos

% 11 – Calcular as componentes das deformações nas direções x e y (ex, ey, Yxy) para cada posição z;
% Sugestão: defina o vetor uma matriz StrainXY = zeros(3,1,zpos) para que ela receba os valores das deformações de forma que:
% StrainXY = zeros(3,1,1) é uma matriz de 3 linhas e uma coluna com os valores das deformações no ponto zpos = 1;

StrainXY = zeros(3, 1, totalPoints);

fprintf(' Results for Strains in XY directions')

for i = 1:n

    for j = 1:npl
        StrainXY(:, 1, j + (i - 1) * npl) = e0 + zpos(j + (i - 1) * npl) * kapa;
        format short e
        fprintf('layer = %d\n', i);
        fprintf('angle = %g\n', orientation(i));
        fprintf('z coordinate = %e\n', zpos(j + (i - 1) * npl));
        fprintf('StrainXY =\n');
        fprintf('  %e %e %e\n', StrainXY(:, 1, j + (i - 1) * npl));

    end

end

% 12 – Calcular as componentes das tensões nas direções x e y (Sx, Sy, Txy) para cada posição z;
% Sugestão: defina o vetor uma matriz StressXY = zeros(3,1,zpos) para que ela receba os valores das deformações de forma que:
% StressXY = zeros(3,1,1) é uma matriz de 3 linhas e uma coluna com os valores das deformações no ponto zpos = 1;

StressXY = zeros(3, 1, totalPoints);

fprintf(' Results for Stresses in XY directions')

for i = 1:n

    for j = 1:npl
        StressXY(:, 1, j + (i - 1) * npl) = Q(:, :, i) * StrainXY(:, 1, j + (i - 1) * npl);
        format short e
        fprintf('layer = %d\n', i);
        fprintf('angle = %g\n', orientation(i));
        fprintf('z coordinate = %e\n', zpos(j + (i - 1) * npl));
        fprintf('StressXY =\n');
        fprintf('  %e %e %e\n', StressXY(:, 1, j + (i - 1) * npl));
    end

end

% 13 – Calcular as componentes das tensões nas direções 1 e 2 (S1, S2, T12) para cada posição z;
% Sugestão: defina o vetor uma matriz Stress12 = zeros(3,1,zpos) para que ela receba os valores das deformações de forma que:
% Stress12 = zeros(3,1,1) é uma matriz de 3 linhas e uma coluna com os valores das deformações no ponto zpos = 1;

Stress12 = zeros(3, 1, totalPoints);

fprintf(' Results for Stresses in 12 directions')

for i = 1:n

    for j = 1:npl
        Stress12(:, 1, j + (i - 1) * npl) = inv(L(:, :, i)') * StressXY(:, 1, j + (i - 1) * npl);
        format short e
        fprintf('layer = %d\n', i);
        fprintf('angle = %g\n', orientation(i));
        fprintf('z coordinate = %e\n', zpos(j + (i - 1) * npl));
        fprintf('Stress12 =\n');
        fprintf('  %e %e %e\n', Stress12(:, 1, j + (i - 1) * npl));
    end

end

S1 = Stress12(1, 1, :)
S2 = Stress12(2, 1, :)
T12 = Stress12(3, 1, :)

% 14 – Calcular as componentes das deformações nas direções 1 e 2 (e1, e2, Y12) para cada posição z;
% Sugestão: defina o vetor uma matriz Strain12 = zeros(3,1,zpos) para que ela receba os valores das deformações de forma que:
% Strain12 = zeros(3,1,1) é uma matriz de 3 linhas e uma coluna com os valores das deformações no ponto zpos = 1;

Strain12 = zeros(3, 1, totalPoints);

fprintf(' Results for Strains in 12 directions')

for i = 1:n

    for j = 1:npl
        Strain12(:, 1, j + (i - 1) * npl) = L(:, :, i) * StrainXY(:, 1, j + (i - 1) * npl);
        format short e
        fprintf('layer = %d\n', i);
        fprintf('angle = %g\n', orientation(i));
        fprintf('z coordinate = %e\n', zpos(j + (i - 1) * npl));
        fprintf('Strain12 =\n');
        fprintf('  %e %e %e\n', Strain12(:, 1, j + (i - 1) * npl));
    end

end

e1 = Strain12(1, 1, :)
e2 = Strain12(2, 1, :)
Y12 = Strain12(3, 1, :)

% 15 – Avaliar para cada posição z o critério de máxima tensão e caso haja falha, indique o modo de falha como foi solicitado na P2;

MT = [Xet, Yet, Xec, Yec, Se];

for i = 1:totalPoints

    if S1(i) > MT(1)
        fprintf(['posicao = %g\n'], zpos(i))
        fprintf(['flag = falha por tração na direção da fibra (%d)\n'], i)
    end

    if S1(i) < -MT(3)
        fprintf(['posicao = %g\n'], zpos(i))
        fprintf(['flag = falha por compressão na direção da fibra (%d)\n'], i)
    end

    if S2(i) > MT(2)
        fprintf(['posicao = %g\n'], zpos(i))
        fprintf(['flag = falha por tração na direção da matriz (%d)\n'], i)
    end

    if S2(i) < -MT(4)
        fprintf(['posicao = %g\n'], zpos(i))
        fprintf(['flag = falha por compressão na direção da matriz (%d)\n'], i)
    end

    if abs(T12(i)) > MT(5)
        fprintf(['posicao = %g\n'], zpos(i))
        fprintf(['flag = falha por cisalhamento no plano (%d)\n'], i)
    end

end

% 16 – Avaliar para cada posição z o critério de máxima deformação e caso haja falha, indique o modo de falha como foi solicitado na P2;

MD = [Xet / E1, Yet / E2, Xec / E1, Yec / E2, Se / G12];

for i = 1:totalPoints

    if e1(i) > MD(1)
        fprintf(['posicao = %g\n'], zpos(i))
        fprintf(['flag = falha por tração na direção da fibra (%d) MAX STRAIN'], i)
    end

    if e1(i) < -MD(3)
        fprintf(['posicao = %g\n'], zpos(i))
        fprintf(['flag = falha por compressão na direção da fibra (%d) MAX STRAIN'], i)
    end

    if e2(i) > MD(2)
        fprintf(['posicao = %g\n'], zpos(i))
        fprintf(['flag = falha por tração na direção da matriz (%d) MAX STRAIN'], i)
    end

    if e2(i) < -MD(4)
        fprintf(['posicao = %g\n'], zpos(i))
        fprintf(['flag = falha por compressão na direção da matriz (%d) MAX STRAIN'], i)
    end

    if abs(Y12(i)) > MD(5)
        fprintf(['posicao = %g\n'], zpos(i))
        fprintf(['flag = falha por cisalhamento no plano (%d)\n'], i)
    end

end

% 17 – Avaliar para cada posição z o critério de Tsai-Hill conforme solicitado na P2;

for i = 1:totalPoints

    if S1(i) > 0
        X = Xet;
    else
        X = Xec;
    end

    if S2(i) > 0
        Y = Yet;
    else
        Y = Yec;
    end

    tsai_hill = (S1(i) / X)^2 - (S1(i) * S2(i)) / (X^2) + (S2(i) / Y)^2 + (T12(i) / Se)^2;

    fprintf(['crit = %g\n'], tsai_hill)
    fprintf(['posicao = %g\n'], zpos(i))

    if tsai_hill < 1

        fprintf(['flag = ok segundo TSAI-HILL\n'])
    else
        fprintf(['flag = falha segundo TSAI-HILL\n'])
    end

end

% 18 – Avaliar para cada posição z o critério de Hoffmann conforme solicitado na P2;

for i = 1:totalPoints

    if S1(i) > 0
        X = Xet;
    else
        X = Xec;
    end

    if S2(i) > 0
        Y = Yet;
    else
        Y = Yec;
    end

    hoffmann = (S1(i)^2) / (Xec * Xet) - (S1(i) * S2(i)) / (Xec * Xet) + (S2(i)^2) / (Yec * Yet) - ((Xet - Xec) / (Xec * Xet)) * S1(i) - ((Yet - Yec) / (Yec * Yet)) * S2(i) + (T12(i)^2) / (Se^2);

    fprintf(['crit = %g\n'], hoffmann)
    fprintf(['posicao = %g\n'], zpos(i))

    if hoffmann < 1
        fprintf(['flag = ok segundo HOFFMANN\n'])
    else
        fprintf(['flag = falha segundo HOFFMANN\n'])
    end

end

% 19 – Avaliar para cada posição z o critério de Tsai-Wu conforme solicitado na P2;

for i = 1:totalPoints

    F1 = 1 / Xet - 1 / Xec;
    F2 = 1 / Yet - 1 / Yec;
    F11 = 1 / (Xet * Xec);
    F22 = 1 / (Yet * Yec);
    F66 = 1 / (Se^2);
    F12 = -sqrt(F11 * F22) / 2;

    tsai_wu = F1 * S1(i) + F2 * S2(i) + F11 * S1(i)^2 + F22 * S2(i)^2 + F66 * T12(i)^2 + 2 * F12 * S1(i) * S2(i);

    fprintf(['crit = %g\n'], tsai_wu)
    fprintf(['posicao = %g\n'], zpos(i))

    if tsai_wu < 1
        fprintf(['flag = ok segundo TSAI-WU com F12 = -(sqrt(F11*F22))/2\n'])
    else
        fprintf(['flag = falha segundo TSAI-WU com F12 = -(sqrt(F11*F22))/2\n'])
    end

end

for i = 1:totalPoints
    F1 = 1 / Xet - 1 / Xec;
    F2 = 1 / Yet - 1 / Yec;
    F11 = 1 / (Xet * Xec);
    F22 = 1 / (Yet * Yec);
    F66 = 1 / (Se^2);
    F12 = 0;

    tsai_wu_2 = F1 * S1(i) + F2 * S2(i) + F11 * S1(i)^2 + F22 * S2(i)^2 + F66 * T12(i)^2 + 2 * F12 * S1(i) * S2(i);

    fprintf(['crit = %g\n'], tsai_wu_2)
    fprintf(['posicao = %g\n'], zpos(i))

    if tsai_wu_2 < 1
        fprintf(['flag = ok segundo TSAI-WU com F12 = 0\n'])
    else
        fprintf(['flag = falha segundo TSAI-WU com F12 = 0\n'])
    end

end
