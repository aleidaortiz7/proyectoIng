clc;
close all;
clear all;
imagenLeer = "mSimetrica.png";

%% PARAMETROS GENERALES (GENERACION DE TODOS LOS CASOS)
generacionAutomatica = 1;
siemprePropagar = 1;

mostrarImagen = 1;
guardarImagen = 0;
mostrarAnimacion = 0;

maxAnchura = 0.1:0.1:0.5;
maxAltura = 0.1:0.1:0.5;
zConjuntoPropagar = (1:5);

% Cantidad de puntos
N = 500;

% Porcentajes de limpieza y curva de decaimiento. Es decir, 
% cuantos porcentaje de pixeles se toman de cada extremo para
% limpiar y para generar una curva de decaimiento.
porcentajeCurvaDecaimiento = 0.025;
    
%% PARAMETROS PARTICULARES (GENERACION DE SOLO UN CASO)
anchuraParticular = [10];
alturaParticular = [10];
zParticular = [10];

if generacionAutomatica == 0
    maxAnchura = anchuraParticular;
    maxAltura = alturaParticular;
    zConjuntoPropagar = zParticular;
end

%% PARAMETROS FASES 1 A 3 Y AJUSTES
for zPropagar = zConjuntoPropagar
    for anchura = maxAnchura
        for altura = maxAltura
            
            % Variable que determina si es un soliton o no
            solitonDescartado = 0;
            
            L = maxAnchura(end);
            
            % Limites inferiores y superiores para los ejes 'x' y 'y'      
            xLimInf = - anchura / 2;
            xLimSup = anchura / 2;
            yLimInf = 0;
            yLimSup = altura;
            
            %% FASE 1 - GENERACION
            % Leer imagen y convertirla a una matriz de escala de grises
            
            % Dicha matriz va de 0 (negro) a blanco (255)
            matPotencialSoliton = imread(imagenLeer);
            
            % En el caso de que la matriz sea a color, tiene 3 dimensiones,
            % por lo que se necesita usar el comando rgb2gray para extrar
            % solo las dimensiones del blanco y negro
            dimColoresImagen = ndims(matPotencialSoliton);
            if dimColoresImagen == 3
              matPotencialSoliton = rgb2gray(matPotencialSoliton); 
            end
            
            % Interpretar la matriz de escalas de grises como un vector de valores 'y'
            
            % Usar el comando flipud para obtener los valores minimos a partir de 
            % la esquina inferior izquierda. De lo contrario, se regresan los min
            % empezando desde la esquina superior izquierda.
            matPotencialSoliton = flipud(matPotencialSoliton);
            
            % Obtener las intensidades minimas y sus posiciones en el eje 'y' para cada
            % columna de la matriz de grises. Filtrar valores que no sean parecidos al negro (0)
            [intensidadesMinimas, yPosMin] = min(matPotencialSoliton);
            yPosMin(intensidadesMinimas >= 10) = 0;
            
            % Encontrar primer y ultima posicion que sean parte del dibujo
            posPrimer = find(yPosMin, 1);
            posUltima = find(yPosMin, 1, 'last');
            
            % Para cada extremo, limpiar un porcentaje de la imagen
            % En el caso del lado izquierdo, se mandan en orden inverso para poder usar
            % la misma funcion sin hacerle muchos cambios. Al final se invierten
            pixelesLimpiar = round((posUltima - posPrimer + 1) * porcentajeCurvaDecaimiento);
            yPosMin =  limpiarExtremos(yPosMin, pixelesLimpiar, posPrimer, posUltima);
            
            % Normalizar el eje 'y' segun los parametros
            yPosMin = rescale(yPosMin, yLimInf, yLimSup);
            
            % Tomar parte significativa de la imagen (no espacio blanco) y extenderla por 3
            
            % Obtener espacio significativo y generar vectores de ello
            tamSignificativo = posUltima - posPrimer;
            
            % Colocar los puntos significativos en el segundo tercio de la imagen
            % 
            xPos = linspace(3 * xLimInf, 3 * xLimSup, tamSignificativo * 3);
            
            yPosMinExt = zeros(1, 3 * tamSignificativo);
            yPosMinExt(tamSignificativo : 2 * tamSignificativo) = yPosMin(posPrimer : posUltima);
            
            % Crear curva de decaimiento para los extremos de la parte significativa
            % Recalcular primer y ultimo valor diferentes a 0
            posPrimer = find(yPosMinExt, 1);
            posUltima = find(yPosMinExt, 1, 'last');
            
            % Tamano de los espacios en ambos lados
            espacioIzq = posPrimer;
            espacioDer = size(yPosMinExt, 2) - posUltima + 1;
            pixelesDecaer = round((posUltima - posPrimer + 1) * porcentajeCurvaDecaimiento);
            
            % Se crean las lineas de decaimiento para cada uno de los lados
            valorDecae = yPosMinExt(posUltima - pixelesDecaer);
            lambda = - log(yPosMinExt(posUltima) / yPosMinExt(posUltima - pixelesDecaer));
            
            curvaDecaimientoDer = valorDecae * exp(- lambda .* (0:(espacioDer - 1)));
            curvaDecaimientoIzq = fliplr(curvaDecaimientoDer);
            
            % Se insertan las curvas en el vector de Y
            yPosMinExt(pixelesDecaer : pixelesDecaer + size(curvaDecaimientoIzq, 2) - 1) = curvaDecaimientoIzq; 
            yPosMinExt((posUltima - pixelesDecaer) : (posUltima - pixelesDecaer) + size(curvaDecaimientoDer, 2) - 1) = curvaDecaimientoDer;  
            
            % Reescalamiento y graficación 
            xPosReescale = linspace(3 * xLimInf, 3 * xLimSup, N);
            
            % Se interpolan los puntos para que la imagen leida tenga los mismo puntos
            % que el codigo elaborado para generar y propagar
            yPosReescale = interp1(xPos, yPosMinExt, xPosReescale);
            yPosReescale(isnan(yPosReescale)) = 0;
            
            % Guardar resultados
            matPotencialSoliton = [xPosReescale; yPosReescale];
            save matPotencialSoliton.mat matPotencialSoliton;
            
            %% FASE 2 - GENERACION
            
            % Leer datos 
            load matPotencialSoliton.mat matPotencialSoliton;
            
            % Calculo de dx y dz que a su vez determina la cantidad de iteraciones
            dx = L / N;
            dz =  dx ^ 2 / 4;         
            numeroIteraciones = zPropagar / dz;
            
            % Se redondea al multiplo de 100 mas cercano para poder hacer
            % bien el step
            numeroIteraciones = ceil(numeroIteraciones / 100) * 100;
            
            % Step de tal forma que siempre haya 100 capas
            stepRevisarPropagacion = numeroIteraciones / 100;
            
            % Transformada de Fourier de la K
            dk = 2 * pi / L;
            k = linspace(-N / 2, N / 2 - 1, N) .* dk;
            kF = (fftshift(k)).^2;
            
            % Funcion del optical field de una imagen
            U = matPotencialSoliton(2, :);
            
            % Calcular la expresion de voltaje  
            lambda = ones(1, N);
            V = lambda - abs(U) .^ 2 - (1 / 2) * ifft((1i)^2 * kF.* fft(U)) ./ (U + 1e-9);
                   
            matCapas = zeros(N, 100);
            numCapa = 1;
            capaActual = zeros(N, 1);
            capaAnterior = zeros(N, 1);  
            
            errorPromedioPerfiles = 0;
            
            for it = 1 : numeroIteraciones    
                % Transformada de Fourier
                % Transformada de Fourier por la exp(-i k^2 deltaZ)
                %Transformada inversa de Fourier
                TF = fft(U);
                TFExp = exp(-1i .* kF .* dz / 4) .* TF;
                TFExpInv = ifft(TFExp);
            
                % Calcular funcion phi y almacenar el resultado
                U = exp(1i .* dz / 2 .* (abs(TFExpInv) .^ 2 + V)) .* TFExpInv;
                
                if it == 1
                    capaAnterior = abs(U);
                end
                
                % Cada ciertas iteraciones, revisar el error entre la iteracion actual
                % y una previa. Si es mucho el error, quiere decir que no es soliton
                if mod(it, stepRevisarPropagacion) == 0
                                      
                    matCapas(:, numCapa) = abs(U);
                    capaActual = matCapas(:, numCapa);
                                   
                    sumCapaAnterior = sum(capaAnterior);
                    sumCapaActual = sum(capaActual);
                    errorPromedioPerfiles = abs(1 - (sumCapaAnterior / sumCapaActual));
                              
                    % En el caso de que no se quiera propagar siempre y que haya mucho
                    % error, se cancela la operacion
                    if  (errorPromedioPerfiles > 0.15)
                        solitonDescartado = 1;
                        
                        if siemprePropagar == 0
                            it = numeroIteraciones;
                        end                   
                    end
                    
                    capaAnterior = capaActual;
                    numCapa = numCapa + 1;
                                                                  
                end % Termina step de revision
                
            end % Termina iteracion
            
            if mostrarImagen == 1              
                if solitonDescartado == 0
                    nombreImg = ['CONFIRMADO_a_', num2str(anchura),'_b_', num2str(altura), num2str(zPropagar), '.png'];
                    colorPlot = 'Blue';
                else
                    nombreImg = ['DESCARTADO_a_', num2str(anchura),'_b_', num2str(altura), num2str(zPropagar), '.png'];
                    colorPlot = 'Red';
                end
                 
                % PLOT SOLITON
                figure(2);  
                plot(xPos, yPosMinExt, 'Color', colorPlot);
                xlim([(-maxAnchura(end) / 2) (maxAnchura(end) / 2)]);
                ylim([0 maxAltura(end)]);
                etiqueta = ['a = ', num2str(anchura), ' b = ', num2str(altura), ' z = ', num2str(zPropagar)];
                title(etiqueta);      
                hold on;
                
                % PLOT MALLA
                             
                if guardarImagen == 1
                    saveas(gcf, nombreImg);  
                end
            end % Termina mostrar imagen
                            
        end % Termina cambio altura
    end % Termina cambio anchura
end % Termina cambio z

% Guardar datos de la variable matCapas en el archivo matDatos.mat
save matSolitonGenerado.mat matCapas;

%% FASE 3 - GRAFICACION

% Graficar cada una de las capas
if mostrarAnimacion == 1
    % Cargar datos
    load matSolitonGenerado.mat matCapas;

    for i = 1 : 100
        plot(xPos, matCapas(:, i));
        xlim([(-maxAnchura(end) / 2) (maxAnchura(end) / 2)]);
        ylim([0 maxAltura(end)]);
        
        title("Propagacion de un soliton");
        pause(0.05);
    end
end

%% FUNCIONES AUXILIARES
% Limpiar los datos para que solo sigan una direccion
function yValue = limpiarExtremos(yValue, cantidadPixeles, posPrimer, posUltima)
    extremoIzq = flip(yValue(posPrimer : posPrimer + cantidadPixeles - 1));
    extremoDer = yValue(posUltima - cantidadPixeles - 1: posUltima);
    
    % Limpia cada extremo
    extremoIzqLimpio = flip(analizarExtremo(extremoIzq));
    extremoDerLimpio = analizarExtremo(extremoDer);
    
    % Reemplaza en los extremos
    yValue(posPrimer : posPrimer + cantidadPixeles - 1) = extremoIzqLimpio;
    yValue(posUltima - cantidadPixeles - 1: posUltima) = extremoDerLimpio;
end

% Funcion que analiza los extremos de una grafica y basado en su tendencia
% ajusta los datos para que sean consistentes.
% El objetivo de esto es que en los extremos o solo haya subida o solo
% bajada para poder extender correctamente.
function datosLimpios = analizarExtremo(datosExtremo)   
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