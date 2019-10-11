clc
close all


% Cargar datos al workspace, se carga la variable matXY
load matDatos.mat matCapas;

% Para encontrar el mayor de las y
yValues(1,:,:) = matCapas(2,:,:);
maxAll = max(yValues,[],'all');

for i = 1 : 100
    % Extraer informacion y graficar
    x = matCapas(1, :, i);
    y = matCapas(2, :, i);
    figure(1);
    plot(x, y);
    ylim([0 maxAll]);
    title("Propagacion soliton");
    pause(0.05);
end