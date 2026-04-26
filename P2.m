clear all
clc

% 1 - Receber as propriedades da lâmina nas direções 1 e 2 (módulos de elasticidade - E1, E2, coeficiente de Poisson - v12, módulo de cisalhamento - G12, resistências da lâmina -Xt,Xc,Yt,Yc,S6 ou Xet,Xec,Yet,Yec,Se);

% E1 = input('Digite o módulo de elasticidade na direção 1 (E1): ');
% E2 = input('Digite o módulo de elasticidade na direção 2 (E2): ');
% v12 = input('Digite o coeficiente de Poisson (v12): ');
% G12 = input('Digite o módulo de cisalhamento (G12): ');
% Xt = input('Digite a resistência à tração na direção 1 (Xt): ');
% Xc = input('Digite a resistência à compressão na direção 1 (Xc): ');
% Yt = input('Digite a resistência à tração na direção 2 (Yt): ');
% Yc = input('Digite a resistência à compressão na direção 2 (Yc): ');
% S6 = input('Digite a resistência ao cisalhamento (S6): ');

E1 = 53780
E2 = 17930
v12 = 0.2500
G12 = 8620
v21 = 0.083349
Xt = 600
Xc = 600
Yt = 50
Yc = 200
S6 = 120

% 2 - Receber o ângulo theta de orientação da lâmina em relação aos eixos x e y;

% theta = input('Digite o ângulo de orientação da lâmina em relação aos eixos x e y (theta): ');
theta = 30

% 3 - Receber o estado de tensão qualquer (Sx, Sy e Txy);

% Sx = input('Digite a tensão na direção x (Sx): ');
% Sy = input('Digite a tensão na direção y (Sy): ');
% Txy = input('Digite a tensão de cisalhamento (Txy): ');

Sx = 500
Sy = -100
Txy = 300

% 4 - Calcular a partir dos dados fornecidos o estado de tensão nas direções 1 e 2 (S1, S2, T12);

c = cosd(theta);
s = sind(theta);

L = [c^2, s^2, c * s;
    s^2, c^2, -c * s;
    -2 * c * s, 2 * c * s, c^2 - s^2]

LT = inv(L')

V12 = LT * [Sx; Sy; Txy]

S1 = V12(1)
S2 = V12(2)
T12 = V12(3)

% 5 - Calcular as componentes das deformações (strains) nas direções 1 e 2 (e1, e2, Y12);


S = [1/E1, -v12/E1, 0;
    -v12/E1, 1/E2, 0;
    0, 0, 1/G12]

E12 = S * [S1; S2; T12]

e1 = E12(1)
e2 = E12(2)
Y12 = E12(3)

% 6 - Calcular as componentes das deformações nas direções x e y (ex, ey, Yxy);

Q = [E1/(1-v12*v21), v12*E2/(1-v12*v21), 0;
    v12*E2/(1-v12*v21), E2/(1-v12*v21), 0;
    0, 0, G12]

QXY =L' * Q * L 

EXY = QXY * [e1; e2; Y12]

ex = EXY(1)
ey = EXY(2)
Yxy = EXY(3)

% 7 - Avaliar se o critério de máxima tensão foi atingido pelo carregamento fornecido;

% 8 – Dizer se falhou ou não (imprimir a informação) e se falhou, dizer qual foi o modo de falha quando o critério foi atingido. Ex.: “Falhou segundo o critério de máxima tensão, o modo de falha foi tração na direção da fibra (1).”

% 9 - Avaliar se o critério de máxima deformação foi atingido pelo carregamento fornecido;

% 10 - Dizer se falhou ou não (imprimir a informação) e se falhou, dizer qual foi o modo de falha quando o critério foi atingido. Ex.: “Falhou segundo o critério de máxima deformação, o modo de falha foi compressão na direção da fibra (1).”

% 11 - Avaliar se o critério de Tsai-Hill foi atingido pelo carregamento fornecido. Dizer se falhou ou não (imprimir a informação) segundo o critério de Tsai-Hill. Ex.: “Não falhou segundo o critério de Tsai-Hill.”

% 12 - Avaliar se o critério de Hoffmann foi atingido pelo carregamento fornecido. Dizer se falhou ou não (imprimir a informação) segundo o critério de Hoffmann. Ex.: “Falhou segundo o critério de Hoffmann.”

% 13 - Avaliar se o critério de Tsai-Wu foi atingido pelo carregamento fornecido. Dizer se falhou ou não (imprimir a informação) segundo o critério de Tsai-Wu. Ex.: “Falhou segundo o critério de Tsai-Wu.”
