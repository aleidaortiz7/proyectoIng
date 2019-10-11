%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Soliton GENERADOR
%
%  Lee la imagen, la convierte a una matriz de blancos y negros y le
%  agrega un decaimiento para poder generar una matriz que metera en el
%  archivo para propagarlo y generar un soliton
%  
% By: Aleida Myriam Machorro Ortiz
% Date: 09/15/2019
% 
% Variables: imag = imagen 
%            imgGris = imagen en grises
%            imgGrisInv = imagen gris invertida
%            minValue = valor del numero minimo de esa columna 
%            yValue = indices de los valores minimos
%            ren = m = numero de filas
%            col = n = numero de columnas
%            posPrimer = primera posicion diferente de cero
%            posUltima = ultima posicion diferente de cero
%            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Leer la imagen y convertirla a escala de grises
clear all
clc 
close all

% Lee y muestra la imagen leida
img = imread('triangulo.png');
figure(1); imshow(img)

% Convierte y muestra la imagen a escala de grises
% Para el caso de imagenes en blanco y negro (generadas con 
% Matlab), se omite el comando rgb2gray que extrae dimensiones
if ndims(img) == 3
  imgGris = rgb2gray(img);  
else
  imgGris = img;
end

figure(2); imshow(imgGris)

%% Convertir la imagen a un vector 

% Calcular los valores minimo de la matriz de grises
% Se invierten los valores de "y" porque se quiere empezar de la esquina
% inferior izquierda, pero min solo puede empezar de la esquina superior
% izquierda
imgGrisInv = flipud(imgGris);
[minValue, yValue] = min(imgGrisInv);

% Hacer cero los valores alejados al negro(0)
yValue(find(minValue > 10)) = 0;

% Cambiar para probar casos
% Sube va con 1, Baja va con -1
direccion = 1;

% Encontrar la primera y ultima posicion diferente de cero en el vec yValue
posPrimer = find(yValue, 1);
posUltima = find(yValue, 1, 'last');

% Para cada extremo, manda los ultimos N pixeles para limpiar
% En el caso del lado izquierdo, se mandan en orden inverso para poder usar
% la misma funcion sin hacerle muchos cambios. Al final se invierten
yValue =  realizaLimpia(yValue, 10, posPrimer, posUltima);

%% Normalizar el eje y

% Calcular el numero filas y de columnas 
[ren, col] = size(imgGris);

% Normalizar el eje "y" [0,1]
yValue = yValue / col;

%% Interpolacion: Imagen a un tercio

% Encontrar la primera y ultima posicion diferente de cero en el vec yValue
posPrimer = find(yValue, 1);
posUltima = find(yValue, 1, 'last');

% Vector que contendra la informacion
y = zeros(1, 3 * (posUltima - posPrimer));
y((posUltima - posPrimer) : 2 * (posUltima - posPrimer)) = yValue(posPrimer : posUltima);

% Vector de las x
xValue = 1 : 3 * (posUltima - posPrimer);


%% Interpolar el vector y

% Se define el vector x
N = 1024;    L = 20;    dx = L / N;
x = linspace(-N / 2, N / 2 - 1, N) .* dx;

%% Crear curva de decaimiento

% Se recalculan el primer y el ultimo valor diferentes a 0
posPrimer = find(y, 1);
posUltima = find(y, 1, 'last');

% Se obtiene el tam de los espacios vacios en ambos lados
espacioIzq = posPrimer;
% Se suma 1 para compensar la resta del rango y que no falten puntos
espacioDer = size(y, 2) - posUltima + 1;

% Se crean las lineas de decaimiento para cada uno de los lados
A = y(posUltima - 5);
lambda = - log(y(posUltima) / y(posUltima - 5));
%lambda = 0.01;
curvaDecaimientoDer = A * exp(- lambda .* (0:(espacioDer - 1)));
curvaDecaimientoIzq = fliplr(curvaDecaimientoDer);

% Se insertan las curvas en el vector de Y
y(5 : (posPrimer + 5)) = curvaDecaimientoIzq; 
y((posUltima - 5) : (end - 5)) = curvaDecaimientoDer;

%% Reescalamiento y graficación
% Reescalar el vector xValue para que vaya de N/2 * dx
xValue = rescale(xValue,-N / 2 * dx, (N / 2 - 1) * dx);

% Se interpolan los puntos para que la imagen leida tenga los mismo puntos
% que el codigo elaborado para generar y propagar
yFinal = interp1(xValue, y, x);
figure(3)
plot(x,yFinal)

%% Guardar imagen

% Gaurdar datos de la variable matImagen en el archivo matDatosImag.mat
matImagen = [x; yFinal];
save matDatosImag.mat matImagen;

%% Limpiar los datos para que solo sigan una direccion
function yValue = realizaLimpia(yValue, cantidadPixeles, posPrimer, posUltima)
    extremoIzq = flip(yValue(posPrimer : posPrimer + cantidadPixeles - 1));
    extremoDer = yValue(posUltima - cantidadPixeles - 1: posUltima);
    
    % Limpia cada extremo
    extremoIzqLimpio = flip(limpiaExtremo(extremoIzq));
    extremoDerLimpio = limpiaExtremo(extremoDer);
    
    % Reemplaza en los extremos
    yValue(posPrimer : posPrimer + cantidadPixeles - 1) = extremoIzqLimpio;
    yValue(posUltima - cantidadPixeles - 1: posUltima) = extremoDerLimpio;
end

% Funcion que analiza los extremos de una grafica y basado en su tendencia
% ajusta los datos para que sean consistentes.
% El objetivo de esto es que en los extremos o solo haya subida o solo
% bajada para poder extender correctamente.
function datosLimpios = limpiaExtremo(datosExtremo)   
   % Determinar la tendencia de esos datos
   tendenciaDatos = obtenerTendencia(datosExtremo);
   
   % Basado en la tendencia, borrar todo a partir del punto donde ya no
   % siga la tendencia. 
   % E.g. Si esta en aumento, 2 5 2 2 7 se convierte en 2 5 0 0 0
  for i = 2 : length(datosExtremo)
    % Caso donde aumenta, queremos que cada elemento sea mayor al anterior
    if tendenciaDatos == 1
        if datosExtremo(i) <= datosExtremo(i - 1)
           datosExtremo(i : end) = 0;
           break;
        end 
    % Caso donde disminuye, queremos que cada elemento sea menor al anterior
    else
        if datosExtremo(i) >= datosExtremo(i - 1)
           datosExtremo(i : end) = 0;
           break;
        end        
    end
  end  
  
  datosLimpios = datosExtremo;
end

% Funcion que dado unos datos, determina si la tendencia de los datos es a
% aumentar o a disminuir. Esto lo hace basado en el vector de diferencias.
function tendenciaDatos = obtenerTendencia(datos)
   difVec = diff(datos);
   
   % Si hay mas positivos, esta aumentando
   if (difVec > 0) > (difVec < 0)
       tendenciaDatos = 1;
   else
       tendenciaDatos = -1;
   end
end