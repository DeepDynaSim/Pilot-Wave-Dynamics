% MATLAB Code for Solving the 2D Time-Dependent Schrödinger Equation

clear;
clc;

% Parameters
NP = 151;
xmax = 1;
xmin = 0;
dx = (xmax - xmin) / (NP-1);
x = linspace(xmin, xmax, NP);
y = linspace(xmin, xmax, NP);
[X, Y] = meshgrid(x, y);

hbar = 1;
m = 1;

% Calculate the maximum allowable time step based on CFL condition
dt_max = (m * dx^2) / hbar;
dt = dt_max * 0.5;
fprintf('Adjusted time step dt = %.5f to satisfy stability.\n', dt);

TimeEnd = 0.003;
nsteps = round(TimeEnd/dt);

% Initial wave packet parameters
x0 = 0.5;
y0 = 0.5;
rho_x = 0.2;
rho_y = 0.2;
kx = 10;
ky = 10;

% Construct the initial wavefunction (Gaussian)
psi = exp(1i * kx * X - ((X - x0).^2) / (2 * rho_x^2)) ...
    .* exp(1i * ky * Y - ((Y - y0).^2) / (2 * rho_y^2));

% Potential (for free particle, V = 0)
V = zeros(NP, NP);

% Precompute some constants to avoid redundant calculations in the loop
hbar2_over_2m = hbar^2 / (2 * m);
dx2 = dx^2;

% Initialize figure for live animation
figure('Renderer', 'painters', 'Position', [100 100 600 500]);
h = surf(X, Y, abs(psi).^2, 'EdgeColor', 'none');
shading interp;
colormap jet;
xlabel('x');
ylabel('y');
zlabel('|\psi(x,y,t)|^2');
axis([xmin xmax xmin xmax 0 max(abs(psi(:)).^2)]);
caxis([0 max(abs(psi(:)).^2)]);
title('Time evolution of the Gaussian wave packet');

% Initialize video writer
video_filename = 'wave_packet_evolution.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 30; % Adjust the frame rate as needed
open(v);

% Time evolution: Implicit scheme
for t = 1:nsteps
    % Kinetic energy operator in x and y (using central differences)
    T_x = (circshift(psi, [0 -1]) - 2 * psi + circshift(psi, [0 1])) / dx2;
    T_y = (circshift(psi, [-1 0]) - 2 * psi + circshift(psi, [1 0])) / dx2;
    
    % Hamiltonian applied to psi
    H_psi = -hbar2_over_2m * (T_x + T_y) + V .* psi;
    
    % Time evolution (Crank-Nicolson Method)
    psi = psi - 1i * dt * H_psi / hbar;
    
    % Update the plot with the new wave function
    if mod(t, nsteps/100) == 0
        set(h, 'ZData', abs(psi).^2);
        title(['Time evolution of the Gaussian wave packet, t = ', num2str(t*dt)]);
        drawnow;

        % Write the current frame to the video
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
end

% Close the video writer
close(v);

% Final visualization: compare the initial and final wave packets
figure;
subplot(1, 2, 1);
surf(X, Y, abs(psi).^2, 'EdgeColor', 'none');
shading interp;
colormap jet;
title('Final wave packet');
xlabel('x');
ylabel('y');
zlabel('|\psi(x,y,t)|^2');
axis([xmin xmax xmin xmax 0 max(abs(psi(:)).^2)]);
caxis([0 max(abs(psi(:)).^2)]);

subplot(1, 2, 2);
initial_psi = exp(1i * kx * X - ((X - x0).^2) / (2 * rho_x^2)) ...
    .* exp(1i * ky * Y - ((Y - y0).^2) / (2 * rho_y^2));
surf(X, Y, abs(initial_psi).^2, 'EdgeColor', 'none');
shading interp;
colormap jet;
title('Initial wave packet');
xlabel('x');
ylabel('y');
zlabel('|\psi(x,y,t)|^2');
axis([xmin xmax xmin xmax 0 max(abs(initial_psi(:)).^2)]);
caxis([0 max(abs(initial_psi(:)).^2)]);