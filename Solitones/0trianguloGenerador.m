% Size and padding of image to be created
sizeImg = 1024;
padding = 250;

matImg = createTriangle(sizeImg, padding);
imshow(matImg);


% Receives sizes, padding and return a matrix representing a triangle
function matImg = createTriangle(sizeImg, padding)
    
    % Multiply by 255 to obtain white background
    matImg = ones(sizeImg, sizeImg) * 255;
   
    % Rows begin in the bottom part of image
    posRow = (sizeImg - padding); 
    
    % Cols begging in left part of image, but need to mirror to paint the
    % other side
    posCol = padding + 1;
    posColMirror = (sizeImg - padding);
    
    % While we don't reach hafl ot image
    while posCol <= (sizeImg / 2)
        % Set pixels to black
        matImg(posRow, posCol) = 0;
        matImg(posRow, posColMirror) = 0;
        
        % Move a row up
        posRow = posRow - 1;

        % Move a column to the left, and mirror that moving the other to
        % the right
        posCol = posCol + 1;
        posColMirror = posColMirror - 1;
    end
end