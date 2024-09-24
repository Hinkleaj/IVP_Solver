% Define the analytical solution
y_exact = @(t) 2*t - 1/2 + 3*exp(-4*t)/2;

% Set up the time range for the plot
t = linspace(0,4,200);

% Compute the analytical solution at the time steps used by the ODE solver
y_exact_vals = y_exact(t);

% Plot the analytical solution
figure
plot(t, y_exact_vals, '-.', 'Color', [0 0.45 0.74], 'LineWidth', 1, 'Marker', '.', 'MarkerSize', 1, 'MarkerEdgeColor', [0.85 0.33 0.1], 'MarkerFaceColor', [0.85 0.33 0.1]);
title('Analytical solution to y′(t) = -4y(t) + 8t, y(0) = 1');
xlabel('t', 'FontSize', 12);
ylabel('y(t)', 'FontSize', 12);
axis([0 4 0 8]);