% MATLAB Code for Simulating 2D TDSE with Gaussian Wave Packet and Cylindrical Potential

% Parameters
NP = 151;                  % Number of grid points
rho = 0.06;                % Initial width of the wave packet
x0 = 0.4;                  % Initial position of the wave packet (x)
y0 = 0.4;                  % Initial position of the wave packet (y)
kx = 120;                  % Wavenumber in x-direction
ky = 120;                  % Wavenumber in y-direction
xmax = 1;                  % Maximum x dimension of the grid
xmin = 0;                  % Minimum x dimension of the grid
ymax = 1;                  % Maximum y dimension of the grid
ymin = 0;                  % Minimum y dimension of the grid
dx = 1/(NP-1);             % Spatial step size (same for x and y)
dt = 0.5 * (dx^2);         % Time step size
lambda = (dx^2) / dt;      % CFL-like parameter

% Potential parameters
ellip = 1;                 % Ellipticity parameter (1 = cylindrical)
potheight = 100000;        % Height of the potential barrier

% Initialize the grid
x = linspace(xmin, xmax, NP);
y = linspace(ymin, ymax, NP);
[X, Y] = meshgrid(x, y);

% Initialize the wave function (Gaussian wave packet)
Inwave = (1/(rho*pi)) * exp(1i*kx*X - ((X-x0).^2)/(2*rho^2)) ...
                  .* exp(1i*ky*Y - ((Y-y0).^2)/(2*rho^2));

% Initialize the potential (cylindrical potential barrier)
V1 = arrayfun(@(x, y) potheight * (norm([(x-0.55) ellip*(y-0.55)]) < 0.05), X, Y);

% Initialize variables
psi = zeros(NP, NP);
phi = Inwave;
chi = zeros(NP, NP);

alpha1 = zeros(NP, NP);
gamma1 = zeros(NP, NP);
alpha2 = zeros(NP, NP);
gamma2 = zeros(NP, NP);
beta1 = zeros(NP, NP);
beta2 = zeros(NP, NP);
phistar = zeros(NP, NP);
psistar = zeros(NP, NP);
eta = zeros(NP, NP);

ap = -1i * lambda;
apc = 1i * lambda;
az = 1 + 2 * 1i * lambda;
azc = 1 - 2 * 1i * lambda;

% Time evolution
TimeEnd = 80;  % Number of time steps
psinew1 = cell(1, TimeEnd);  % Store the wave function at each time step

for dti = 1:TimeEnd
    % Perform calculations as per the discretized TDSE
    % First pass (alpha1 and gamma1 updates)
    for m = 1:NP
        gamma1(NP-1, m) = -1/az;
        for j = NP-1:-1:2
            alpha1(j-1, m) = gamma1(j, m) * ap;
            gamma1(j-1, m) = -1 / (az + ap * alpha1(j-1, m));
        end
    end

    % Second pass (alpha2 and gamma2 updates)
    for j = 1:NP
        gamma2(j, NP-1) = gamma1(NP-1, j);
        for m = NP-1:-1:2
            alpha2(j, m-1) = gamma2(j, m) * ap;
            gamma2(j, m-1) = -1 / (az + ap * alpha2(j, m-1));
        end
    end

    % Forward passes for phi and chi
    for j = 2:NP-1
        for m = 2:NP-1
            chi(j, m) = apc * (phi(j, m+1) + phi(j, m-1)) + azc * phi(j, m);
        end
    end

    for m = 2:NP-1
        beta1(NP-1, m) = 0;
        for j = NP-1:-1:2
            beta1(j-1, m) = gamma1(j, m) * (ap * beta1(j, m) - chi(j, m));
        end
    end

    for m = 2:NP-1
        phistar(1, m) = 0;
        phistar(NP, m) = 0;
        for j = 1:NP-1
            phistar(j+1, m) = alpha1(j, m) * phistar(j, m) + beta1(j, m);
        end
    end

    % Apply the potential and update psistar
    for j = 1:NP
        for m = 2:NP-1
            psistar(j, m) = exp(-0.5 * 1i * dt * V1(j, m)) * phistar(j, m);
        end
    end

    % Forward passes for psi and eta
    for j = 2:NP-1
        for m = 2:NP-1
            eta(j, m) = apc * (psistar(j+1, m) + psistar(j-1, m)) + azc * psistar(j, m);
        end
    end

    for j = 2:NP-1
        beta2(j, NP-1) = 0;
        for m = NP-1:-1:2
            beta2(j, m-1) = gamma2(j, m) * (ap * beta2(j, m) - eta(j, m));
        end
    end

    for j = 2:NP-1
        psi(j, 1) = 0;
        for m = 1:NP-1
            psi(j, m+1) = alpha2(j, m) * psi(j, m) + beta2(j, m);
        end
    end

    % Apply the potential and update phi
    for m = 2:NP-1
        for j = 2:NP-1
            phi(j, m) = exp(-0.5 * 1i * dt * V1(j, m)) * psi(j, m);
        end
    end

    psinew1{dti} = phi;  % Store the updated wave function
end

% Plot the results similar to Mathematica code
figure;
for i = 1:4
    subplot(2, 2, i);
    time_index = round(linspace(1, TimeEnd, 4));
    imagesc(abs(psinew1{time_index(i)}).^2 + 2 * V1 / potheight);
    colorbar;
    axis equal tight;
    xlabel('x');
    ylabel('y');
    title(['Time Step: ', num2str(time_index(i))]);
end

sgtitle('Wave Function Interaction with Cylindrical Potential');

% Optionally, create an animation
figure;
filename = '2DCylindricalPotential3D.gif';
for i = 1:TimeEnd
    imagesc(abs(psinew1{i}).^2 + 2 * V1 / potheight);
    colorbar;
    axis equal tight;
    title(['Time Step: ', num2str(i)]);
    drawnow;

    % Capture the plot as an image
    frame = getframe(gcf);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);

    % Write to the GIF File
    if i == 1
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

disp('Animation saved as 2DCylindricalPotential3D.gif');
