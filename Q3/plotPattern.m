function plotPattern(pattern, theta, phi)
    % plotPattern: Visualizes a 3D gain pattern in Cartesian coordinates.
    %
    % Inputs:
    %   pattern - Matrix of gain values (linear scale)
    %   theta - Vector of theta values (degrees)
    %   phi - Vector of phi values (degrees)

    % Convert spherical coordinates to Cartesian coordinates
    [phi_grid, theta_grid] = meshgrid(deg2rad(phi), deg2rad(theta)); % Convert to radians
    r = pattern; % Assume pattern is the radius (gain in linear scale)
    x = r .* sin(theta_grid) .* cos(phi_grid);
    y = r .* sin(theta_grid) .* sin(phi_grid);
    z = r .* cos(theta_grid);

    % Create 3D surface plot
    figure;
    surf(x, y, z, 10 * log10(r), "EdgeAlpha", .4); % Use r in dB for coloring
    colormap(.7*ones(1,3)); % Grayscale colormap

    % Add labels and title
    xlabel('X (Normalized)');
    ylabel('Y (Normalized)');
    zlabel('Z (Normalized)');
    title('3D Gain Pattern in Cartesian Coordinates (Grayscale)');

    % Set axis properties
    axis equal; % Equal scaling for all axes
    grid on;
    % view(3); % 3D perspective
end
