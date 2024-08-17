% MATLAB Code for Visualizing Square Well, Elliptical, and Cylindrical Potential Barriers

% Parameters
NP = 151;               % Number of grid points
potHeight = 10000;      % Height of the potential barrier

% Define the grid
x = linspace(0, 1, NP);
y = linspace(0, 1, NP);
[X, Y] = meshgrid(x, y);

% Square Well Potential Barrier
V1 = zeros(NP, NP);
V1(X > 0.5 & X < 0.7 & Y > 0.2 & Y < 0.8) = potHeight;

% Elliptical Potential Barrier
V2 = zeros(NP, NP);
V2(sqrt((X - 0.5).^2 + (0.5 * (Y - 0.5)).^2) < 0.1) = potHeight;

% Cylindrical (Circular) Potential Barrier
V3 = zeros(NP, NP);
V3(sqrt((X - 0.5).^2 + (Y - 0.5).^2) < 0.1) = potHeight;

% Visualize the potentials
figure;
subplot(1,3,1);
imagesc(x, y, V1);
axis xy;
colorbar;
title('Square Well Potential');
xlabel('x');
ylabel('y');

subplot(1,3,2);
imagesc(x, y, V2);
axis xy;
colorbar;
title('Elliptical Potential');
xlabel('x');
ylabel('y');

subplot(1,3,3);
imagesc(x, y, V3);
axis xy;
colorbar;
title('Cylindrical (Circular) Potential');
xlabel('x');
ylabel('y');

% Adjust the figure properties
sgtitle('Comparison of Different Potential Barriers');
set(gcf, 'Position', [100, 100, 1200, 400]);  % Adjust figure size
