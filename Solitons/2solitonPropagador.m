% Con una z de 20 y con 20972 iteraciones
% Generar y propagar un perfil de onda con un voltaje 
clear all
close all

N = 1024;
L = 20;
dx = L / N;
dz =  dx ^ 2 / 4;
x = linspace(-N / 2, N / 2 - 1, N) .* dx;
dk = 2 * pi / L;
k = linspace(-N / 2, N / 2 - 1, N) .* dk;
% Transformada de Fourier de la K
kF = (fftshift(k)).^2;

% Pregunta cuanto es lo que se propagara el haz en Z
z = input('Longitud a propagar en Z: ')
% Numero de iteraciones
noDeIt = z / dz;

% Funcion del optical field de una imagen
load matDatosImag.mat matImagen;
U = matImagen(2, :);
% indUsados = matImagen(1,:);
% x = x(indUsados);
% k = k(indUsados);

%U = 2 * sech(x);
%U = exp(- x .^ 2);
%plot(x,U)

% Calcular la expresion de voltaje  
lambda = ones(1,N);
V = lambda - abs(U) .^ 2 - (1 / 2) * ifft((1i)^2 * kF.* fft(U)) ./ (U + 1e-9);
%figure(2)
%plot(x,U)

% Esto es para que no s2ea un soliton/Probar que funciona mi codigo
%U = U * 2 ;

% Se calcula el step para solo guardar 100 capas
step = floor(noDeIt / 100);

% Numero de capas guardadas
noCapas = 0;
% Matriz donde se guardara las capas
matCapas = zeros(2,N,100);
% Numero de comparaciones para ir checando que le perfil se mantenga
numComparaciones = 4;
% Step entre el numero de capas para checar que el perfil se mantenga
stepPerfilEstable = 100 / numComparaciones;
% Contador de la posicion de las compraciones
contComparacionesPerfil = stepPerfilEstable;
% Calcula el error entre los perfiles
errorPerfiles = 0;

for r = 1 : noDeIt
    % Transformada de Fourier
    TF = fft(U);

    % Transformada de Fourier por la exp(-i k^2 deltaZ)
    TFExp = exp(-1i .* kF .* dz / 4) .* TF;
   
    %Transformada inversa de Fourier 
    TFExpInv = ifft(TFExp);

    % La funcion phi
    U = exp(1i .* dz / 2 .* (abs(TFExpInv) .^ 2 + V)) .* TFExpInv;
    
    if mod(r,step) == 0
        noCapas = noCapas + 1
        matCapas(1,:,noCapas) = x; 
        matCapas(2,:,noCapas) = abs(U); 
        if noCapas == contComparacionesPerfil
            errorPerfiles = sum(abs(matCapas(2,:,contComparacionesPerfil) - matCapas(2,:,contComparacionesPerfil - stepPerfilEstable + 1))) / N;
            contComparacionesPerfil = contComparacionesPerfil + stepPerfilEstable;
            % Si el error es mas grande que 1x10^-3, se para el codigo
            if errorPerfiles > 1e-3
                disp('NO ES UN SOLITON')
                break
            end
        end
    end
    
end

% Guardar datos de la variable matCapas en el archivo matDatos.mat
save matDatos.mat matCapas;