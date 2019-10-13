clc;
close all;
clear all;

% Define folder to analyze and name of text file
folderName = "IMAGE-DMD-A1-0-1.A2-0.5-AYA";
%folderName = "IMAGE-DMD-A1-0-1.A2-1-AYA";

% Move to folder and get the list of all bmp files
folderFiles = dir('*.bmp');
countFiles = size(folderFiles, 1);

% Output variables to store mean intensity of each image and each radius
outputIntensity = zeros(2, countFiles);
outputDifference = zeros(1, countFiles);

% / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /

fileToAnalyze = folderFiles(1);
[matLight, imageLight, imageSizeX, imageSizeY] = convertImage(fileToAnalyze);

% Max radius to test mean intensity 
maxRadiusToTest = 10;

% Variables to store max difference, radius that created it and the circles
% in the screen that correspond to that area.
maxDifference = 0;
bestRadius = 0;
bestCircles = zeros(imageSizeX, imageSizeY);

% / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /

figure(1);
imshow(imageLight);

% Select 2 points to analyze, its coordinates will be stored
[xCoord, yCoord] = ginput(2);
yCoord = round(yCoord); 
xCoord = round(xCoord); 

% For each of the radiuses to test, calculate the intensity in that region
for radiusTest = 1 : maxRadiusToTest
    [meanIntensity, circleArea, intDifference] = calculateIntensities(matLight, xCoord, yCoord, radiusTest, imageSizeX, imageSizeY);

    % Check if this difference is best to what we had before
    if intDifference > maxDifference
        maxDifference = intDifference;
        bestRadius = radiusTest;
        bestCircles = circleArea;
    end
end

% Store that in the output matrices, meanIntensity is transposed because we
% want to have a set of data per column
outputIntensity(:, 1) = meanIntensity';
outputDifference(1) = intDifference;

% Comment or uncomment to display resultss
displayResults(intDifference, bestRadius, xCoord, yCoord, bestCircles);

% For each of the reamaining files in the folder, convert the image to an
% intensity matrix, analyze the intensity difference and store those values
for fileIndex = 2 : countFiles
    fileToAnalyze = folderFiles(fileIndex);
    [matLight, imageLight, imageSizeX, imageSizeY] = convertImage(fileToAnalyze);
    
    [meanIntensity, circleAre, intDifference] = calculateIntensities(matLight, xCoord, yCoord, bestRadius, imageSizeX, imageSizeY);
    
    outputIntensity(:, fileIndex) = meanIntensity';
    outputDifference(fileIndex) = intDifference;
end

% Save output values in .mat file
save(folderName + '.mat', 'outputIntensity', 'outputDifference', 'bestRadius', 'xCoord', 'yCoord');

% Calculate gamma 
for i = 1 : countFiles
	gamma(i) = abs(outputIntensity(1,i) - outputIntensity(2,i)) / (outputIntensity(1,i) + outputIntensity(2,i));
end

figure()
xX = 0 : (1 / (countFiles - 1)) : 1;
plot(xX, gamma)

% / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /

% Given a file, open it in the current folder and convert it to an
% intensity matrix
function[matLight, imageLight, imageSizeX, imageSizeY] = convertImage(fileToAnalyze)
    % Load image and obtain matrix of pixel intensity from 0 (black) to 255 (white)
    imageLight = imread(fileToAnalyze.name);
    matLight = rgb2gray(imageLight); 
    [imageSizeY, imageSizeX] = size(matLight);
end

% Given an intensity matrix, a radius and coordinates, analyze the
% intensity in the circle that is formed with that radius, and calculate
% the difference in intensities.
function[meanIntensity, circleArea, intDifference] = calculateIntensities(matLight, xCoord, yCoord, radius, imageSizeX, imageSizeY)
    % Mean intensities and area where these are being obtained
    meanIntensity = zeros(1, 2);
    circleArea = zeros(imageSizeY, imageSizeX);
   
    for i = 1 : 2
        % Create a radius of pixels around the clicked points. 
        % A squared grid is created and then we filter out points that are not
        % inside the area of the circle.
        [colsInImage, rowsInImage] = meshgrid(1 : imageSizeX, 1 : imageSizeY);

        circlePixels = (rowsInImage - yCoord(i)).^2 ... 
                     + (colsInImage - xCoord(i)).^2 ... 
                    <= ((radius.^2)); 

        matIntensity = matLight(circlePixels);
        meanIntensity(i) = mean(matIntensity,'all') / 255;
        
        circleArea(circlePixels) = 255;
    end 
    
    intDifference = abs(meanIntensity(2) - meanIntensity(1));
end

% / / / / / / / / / / / / / / / / / / / / / / / / / / / / / /

function displayResults(intDifference, bestRadius, xCoord, yCoord, bestCircles)    
    % Display the max difference obtained and radius of it
    txtDiff = ['Max diff: ' num2str(intDifference)];
    txtRad = ['Radius: ', num2str(bestRadius)];
    disp(txtDiff);
    disp(txtRad);

    % Also display coordinates of clicked points
    txtPos1 = ['X1: ' num2str(xCoord(1)), ', Y1: ', num2str(yCoord(1))];
    txtPos2 = ['X2: ' num2str(xCoord(2)), ', Y2: ', num2str(yCoord(2))];
    disp(txtPos1);
    disp(txtPos2);
    
    posMarker = [xCoord(1) yCoord(1); xCoord(2) yCoord(2)];
    bestCirclesShow = insertMarker(bestCircles, posMarker,'x','color', 'red','size', 10);
    
    figure(2);
    imshow(bestCirclesShow);
    axis('on', 'image');
end