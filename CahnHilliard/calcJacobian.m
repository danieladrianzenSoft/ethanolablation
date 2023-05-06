function jacobian = calcJacobian(matrix)
% Define the convolution kernels for gradient calculation
dxKernel = [-1 0 1; -2 0 2; -1 0 1];
dyKernel = dxKernel';

% Pad the matrix with edge wrapping
paddedMatrix = padarray(matrix, [1 1], 'circular');

% Perform the convolutions using imfilter function
dx = imfilter(paddedMatrix, dxKernel, 'conv', 'replicate');
dy = imfilter(paddedMatrix, dyKernel, 'conv', 'replicate');

% Crop the convolved matrices back to the original size of the input
dx = imcrop(dx, [2 2 size(matrix, 2) - 1 size(matrix, 1) - 1]);
dy = imcrop(dy, [2 2 size(matrix, 2) - 1 size(matrix, 1) - 1]);

% Calculate the Jacobian by combining the gradient matrices
jacobian = cat(3, dx, dy);

end