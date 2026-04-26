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

S = [1 / E1, -v12 / E1, 0;
    -v12 / E1, 1 / E2, 0;
    0, 0, 1 / G12]

Q = [E1 / (1 - v12 * v21), v12 * E2 / (1 - v12 * v21), 0;
    v12 * E2 / (1 - v12 * v21), E2 / (1 - v12 * v21), 0;
    0, 0, G12]

% 2 - Receber o ângulo theta de orientação da lâmina em relação aos eixos x e y;

% theta = input('Digite o ângulo de orientação da lâmina em relação aos eixos x e y (theta): ');
theta = 30

s = sind(theta)
c = cosd(theta)

L = [c^2, s^2, c * s;
    s^2, c^2, -c * s;
    -2 * c * s, 2 * c * s, c^2 - s^2]                                                                                                     

LT = L'
LTi = inv(LT)

% 3 - Receber o estado de tensão qualquer (Sx, Sy e Txy);

% Sx = input('Digite a tensão na direção x (Sx): ');
% Sy = input('Digite a tensão na direção y (Sy): ');
% Txy = input('Digite a tensão de cisalhamento (Txy): ');

Sx = 500
Sy = -100
Txy = 300

Vxy = [Sx; Sy; Txy]

% 4 - Calcular a partir dos dados fornecidos o estado de tensão nas direções 1 e 2 (S1, S2, T12);

V12 = LTi * Vxy

% 5 - Calcular as componentes das deformações (strains) nas direções 1 e 2 (e1, e2, Y12);

E12 = S * V12

e1 = E12(1);
e2 = E12(2);
Y12 = E12(3);

% 6 - Calcular as componentes das deformações nas direções x e y (ex, ey, Yxy);

QXY = L' * Q * L

EXY = QXY * [e1; e2; Y12]

S1 = V12(1)
S2 = V12(2)
T12 = V12(3)

% Não entendi, ele pede essa resposta no enunciado mas não quer na resposta oficial que ele deu?
ex = EXY(1)
ey = EXY(2)
Yxy = EXY(3)

% 7 - Avaliar se o critério de máxima tensão foi atingido pelo carregamento fornecido;

MT = [Xt, Yt, Xc, Yc, S6]

% 8 – Dizer se falhou ou não (imprimir a informação) e se falhou, dizer qual foi o modo de falha quando o critério foi atingido. Ex.: “Falhou segundo o critério de máxima tensão, o modo de falha foi tração na direção da fibra (1).”

direcao = {"fibra", "matriz"}

for i = 1:2

    if V12(i) < MT(i) && V12(i) > -MT(i + 2)
        % fprintf('Não falhou segundo o critério de máxima tensão na direção (%d).\n', i);
    else

        if V12(i) > MT(i)
            fprintf('Falhou segundo o critério de máxima tensão, o modo de falha foi tração na direção da %s (%d).\n', direcao{i}, i);
        end

        if V12(i) < -MT(i + 2)
            fprintf('Falhou segundo o critério de máxima tensão, o modo de falha foi compressão na direção da %s (%d).\n', direcao{i}, i);
        end

    end

end

if abs(V12(3)) < MT(5)
    % fprintf('Não falhou segundo o critério de máxima tensão no cisalhamento.\n');
else
    fprintf('Falhou segundo o critério de máxima tensão, o modo de falha foi cisalhamento.\n');
end

% 9 - Avaliar se o critério de máxima deformação foi atingido pelo carregamento fornecido;

MD = [Xt / E1, Yt / E2, Xc / E1, Yc / E2, S6 / G12]

% 10 - Dizer se falhou ou não (imprimir a informação) e se falhou, dizer qual foi o modo de falha quando o critério foi atingido. Ex.: “Falhou segundo o critério de máxima deformação, o modo de falha foi compressão na direção da fibra (1).”

for i = 1:2

    if E12(i) < MD(i) && E12(i) > -MD(i + 2)
        % fprintf('Não falhou segundo o critério de máxima deformação na direção (%d).\n', i);
    else

        if E12(i) > MD(i)
            fprintf('Falhou segundo o critério de máxima deformação, o modo de falha foi tração a direção da %s (%d).\n', direcao{i}, i);
        end

        if E12(i) < -MD(i + 2)
            fprintf('Falhou segundo o critério de máxima deformação, o modo de falha foi compressão a direção da %s (%d).\n', direcao{i}, i);

        end

    end

end

% 11 - Avaliar se o critério de Tsai-Hill foi atingido pelo carregamento fornecido. Dizer se falhou ou não (imprimir a informação) segundo o critério de Tsai-Hill. Ex.: “Não falhou segundo o critério de Tsai-Hill.”

if S1 > 0
    X = Xt
else
    X = Xc
end

if S2 > 0
    Y = Yt
else
    Y = Yc
end

tsai_hill = (S1 / X)^2 - (S1 * S2) / (X^2) + (S2 / Y)^2 + (T12 / S6)^2

if tsai_hill > 1
    fprintf('Falhou segundo o critério de Tsai-Hill.\n');
end

% 12 - Avaliar se o critério de Hoffmann foi atingido pelo carregamento fornecido. Dizer se falhou ou não (imprimir a informação) segundo o critério de Hoffmann. Ex.: “Falhou segundo o critério de Hoffmann.”

hoffmann = (S1^2) / (Xc * Xt) - (S1 * S2) / (Xc * Xt) + (S2^2) / (Yc * Yt) - ((Xt - Xc) / (Xc * Xt)) * S1 - ((Yt - Yc) / (Yc * Yt)) * S2 + (T12^2) / (S6^2)

if hoffmann > 1
    fprintf('Falhou segundo o critério de Hoffmann.\n');
end

% 13 - Avaliar se o critério de Tsai-Wu foi atingido pelo carregamento fornecido. Dizer se falhou ou não (imprimir a informação) segundo o critério de Tsai-Wu. Ex.: “Falhou segundo o critério de Tsai-Wu.”

F1 = 1 / Xt - 1 / Xc;
F2 = 1 / Yt - 1 / Yc;
F11 = 1 / (Xt * Xc);
F22 = 1 / (Yt * Yc);
F66 = 1 / (S6^2);
F12 = -sqrt(F11 * F22)/2;

tsai_wu_1 = F1 * S1 + F2 * S2 + F11 * S1^2 + F22 * S2^2 + F66 * T12^2 + 2 * F12 * S1 * S2

if tsai_wu_1 > 1
    fprintf('Falhou segundo Tsai-Wu com F12 = -(sqrt(F11*F22))/2.\n');
end

F12 = 0;

tsai_wu_2 = F1 * S1 + F2 * S2 + F11 * S1^2 + F22 * S2^2 + F66 * T12^2 + 2 * F12 * S1 * S2

if tsai_wu_2 > 1
    fprintf('Falhou segundo Tsai-Wu com F12 = 0.\n');
end