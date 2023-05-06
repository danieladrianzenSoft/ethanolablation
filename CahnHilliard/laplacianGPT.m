function laplacian = laplacianGPT(matrix,dx)
% Define the convolution kernel for Laplacian calculation
kernel = (1/(dx^2))*[0 1 0; 1 -4 1; 0 1 0];

% Pad the matrix with edge wrapping
paddedMatrix = padarray(matrix, [1 1], 'circular');

% Perform the convolution using imfilter function
laplacian = imfilter(paddedMatrix, kernel, 'conv', 'replicate');

% Crop the convolved matrix back to the original size of the input
laplacian = imcrop(laplacian, [2 2 size(matrix, 2) - 1 size(matrix, 1) - 1]);

end