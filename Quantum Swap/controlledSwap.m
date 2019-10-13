%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controlled Swap for qutrits
% Elaborated by: Aleida Machorro Ortiz
% Date: 2/26/2019
% This program calculates Kolmogorov constant for LG with opposite
% topological charge from 1 to M and amplitudes that varies from exp(i0), 
% exp(i2pi/3) and exp(i4pi/3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%
% Close all and clear the variables
close all;
clear all;

% Input values
punt = 1024;  
w0 = 0.5e-3;

% Create vectors in x and y
xmax = 1.5e-3;  ymax = xmax; pm = punt / 2;
x = linspace(-xmax, xmax * (punt - 2) / punt, punt); 
y = linspace(ymax, ymax * (2 - punt) / punt, punt);

[xx, yy] = meshgrid(x, y);
dx = abs(((xmax * (punt - 2) / punt)- (-xmax)) / (punt - 1));
dy = abs(((ymax * (2 - punt) / punt)- (ymax)) / (punt - 1));

% Coordinate transformation
rr = sqrt(xx .^ 2 + yy .^ 2);  phi = atan2(yy, xx);
[nRow, nCol] = size(rr);

% Laguerre z = 0, order n = 0 and order azimutal m;
LG = @(m) rr .* exp(-rr .^ 2 / w0 ^ 2) .* exp(1i * m * phi);

% How much LG do you want?
prompt = 'How much LG do you want? ';
M = input(prompt);

% Genera la tabla de Qutrits
tabQu = GeneraLaTabQu(M);

% Sum f(m,-m) and f(-m,m)
for numMatriz = 1 : 3 ^ (2 * M)
    sumMas = zeros(nRow, nCol);
    sumMenos = zeros(nRow, nCol);
    intVer = zeros(nRow, nCol);
    intHor = zeros(nRow, nCol);
    kAnalitica = zeros(1, 1);

    % Cont columns
    col = 1;
    etaAnalitica = 0;

    for j = 1 : M
        % f(m,-m)
        sumMas(:, :) = sumMas(:, :) + tabQu(numMatriz, col) ...
                          * LG(j) + tabQu(numMatriz, col + 1) * LG(-j);

        % f(-m,m)
        sumMenos(:, :) = sumMenos(:, :) + tabQu(numMatriz, col) ...
                          * LG(-j) + tabQu(numMatriz, col + 1) * LG(j);

        % eta Analitica
        etaAnalitica =  etaAnalitica + abs(tabQu(numMatriz, col) ...
                        - tabQu(numMatriz, col + 1)) ^ 2;

        col = col + 2;
    end

    tabQu(numMatriz, 4 * M + 1) = etaAnalitica / (4 * M);

    % Intensity vertical = 1/2 * (f(m,-m) - f(-m,m))
    intVer = (1 / 2) .* (sumMas - sumMenos) .* (1 / 2) .* conj((sumMas - sumMenos));
    % Intensity horizontal = 1/2 * (f(m,-m) + f(-m,m))
    intHor = (1 / 2) .* (sumMas + sumMenos) .* (1 / 2) .* conj((sumMas + sumMenos));

    % Integral de la potencia
    Pv = sum(sum(intVer(:, :))) * dx * dy;
    Ph = sum(sum(intHor(:, :))) * dx * dy;

    %Ph y Pv normalizadas
    PhNor = Ph ./ (Ph + Pv);
    PvNor = Pv ./ (Ph + Pv);

    % Constante eta
    eta = PvNor;
    tabQu(numMatriz, 4 * M + 2) = eta;
end

nombreArchivo = ['datosM', num2str(M), '.mat'];
save(nombreArchivo, 'tabQu');

% Encontrar el numero de coincidencias
vecValCoinc = 0 : (0.75 / M) : 0.75;
% Se obtiene el pocentaje de coincidencias
porcentajeCoinc = (0 : 1 / M : 1) .* 100;
plot(porcentajeCoinc, vecValCoinc,'o-')
xlabel('% porcentaje de error')
ylabel('Valor de eta^2')

%%%%%%%%%%%% FUNCIONES %%%%%%%%%%%%
% Genera la tabla de qutrits 
function tabQu = GeneraLaTabQu(M) 
    tabQu = zeros(3 ^ (2 * M),(2 * M));
    %Contador de cuantas veces se ha ejecutado el ciclo de la col
    contVecesCol = 0;

    %GENERA LA TABLA DE QUTRITS QUE TENDRAN LAS DIFERENTES INTENSIDADES
    %Me ayuda moverme columna por columna
    for col = (2 * M) : -1 : 1
        %Cuanta cuantas veces se ha ejecutado el ciclo del renglon
        contVecesRenglon = 1;
        %Cuenta en que renglon esta
        contR = 1; 
        %Contador de repeticiones
        contRep = 3 ^ contVecesCol;

        %Mueve renglon por renglon
        for ren =  1 : contRep : 3 ^ (2 * M)
            %Le asigna el valor
            if mod(contVecesRenglon,3) == 1
                %Valor 1 de a
                valEtiqueta = 1;
                val = exp(1i * 0);
            elseif mod(contVecesRenglon, 3) == 2
                %Valor 2 de b
                valEtiqueta = 2;
                val = exp(1i * 2 * pi / 3);
            else
                %Valor 3 de c
                valEtiqueta = 3;
                val = exp(1i * 4 * pi / 3);
            end
            %Genera la repeticiones del valor por columna
            for rep = 1 : contRep
                tabQu(contR, col) = val;
                tabQu(contR, col + 2 * M) = valEtiqueta;
                contR = contR + 1;
            end
            contVecesRenglon = contVecesRenglon + 1;
        end
        contVecesCol = contVecesCol + 1;
    end
end