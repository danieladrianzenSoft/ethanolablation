% Define the initial conditions
nx = 100; % number of points in the x direction
ny = 100; % number of points in the y direction
D = 1e-4; % diffusion coefficient
t_end = 100; % end time
dt = 0.5; % time step
dx = 0.05; % spatial step
dy = 0.05; % spatial step
C = zeros(nx,ny); % initial concentration

% Define the boundary conditions
C(1:5,:) = 1; % left boundary
C(nx,:) = 0; % right boundary
C(:,1) = 0; % bottom boundary
C(:,ny) = 0; % top boundary

% Implement the finite difference method
for t = 0:dt:t_end
    for i = 2:nx-1
        for j = 2:ny-1
            C(i,j) = C(i,j) + D*dt*(C(i+1,j)-2*C(i,j)+C(i-1,j))/(dx^2) ...
                     + D*dt*(C(i,j+1)-2*C(i,j)+C(i,j-1))/(dy^2);
        end
    end
    
    % Plot the concentration at each time step
    imagesc(C);
    colorbar;
    title(sprintf('Time: %0.2f seconds',t));
    xlabel('X');
    ylabel('Y');
    drawnow;
end