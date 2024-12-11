function projectile_gui
    % Create the main figure window
    fig = figure('Name', 'Projectile Motion Solver', ...
                 'NumberTitle', 'off', ...
                 'Position', [300, 300, 600, 500]);

    % Create UI components
    uicontrol('Style', 'text', 'Position', [50, 440, 500, 20], ...
              'String', 'Projectile Motion Solver with Visualization', ...
              'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    uicontrol('Style', 'text', 'Position', [50, 400, 150, 20], ...
              'String', 'Initial Velocity (m/s):', 'HorizontalAlignment', 'left');
    velocityInput = uicontrol('Style', 'edit', 'Position', [200, 400, 150, 25]);

    uicontrol('Style', 'text', 'Position', [50, 360, 150, 20], ...
              'String', 'Distance (m):', 'HorizontalAlignment', 'left');
    distanceInput = uicontrol('Style', 'edit', 'Position', [200, 360, 150, 25]);

    uicontrol('Style', 'text', 'Position', [50, 320, 150, 20], ...
              'String', 'Height (m):', 'HorizontalAlignment', 'left');
    heightInput = uicontrol('Style', 'edit', 'Position', [200, 320, 150, 25]);

    uicontrol('Style', 'text', 'Position', [50, 280, 150, 20], ...
              'String', 'Gravity (m/s^2):', 'HorizontalAlignment', 'left');
    gravityInput = uicontrol('Style', 'edit', 'Position', [200, 280, 150, 25], 'String', '9.8');

    % Results
    resultText = uicontrol('Style', 'text', 'Position', [50, 240, 500, 40], ...
                           'String', '', 'HorizontalAlignment', 'center', ...
                           'FontSize', 12, 'ForegroundColor', 'blue');

    % Axes for the plot
    ax = axes('Parent', fig, 'Position', [0.1, 0.1, 0.8, 0.3]);

    % Compute Button
    uicontrol('Style', 'pushbutton', 'Position', [250, 200, 100, 40], ...
              'String', 'Compute', ...
              'Callback', @computeProjectile);

    % Callback function for the Compute button
    function computeProjectile(~, ~)
        % Read inputs
        v0 = str2double(get(velocityInput, 'String'));
        R = str2double(get(distanceInput, 'String'));
        H = str2double(get(heightInput, 'String'));
        g = str2double(get(gravityInput, 'String'));

        % Validate inputs
        if isnan(v0) || isnan(R) || isnan(H) || isnan(g)
            set(resultText, 'String', 'Invalid input. Please enter valid numbers.');
            return;
        end

        % Solve for launch angle and time using Newton-Raphson
        try
            [theta, time] = projectile_newton_raphson(v0, R, H, g);
            resultString = sprintf('Launch Angle: %.2f°\nTime of Flight: %.2f s', rad2deg(theta), time);
            % Plot trajectory
            plotTrajectory(v0, theta, time, g, ax);
        catch err
            resultString = sprintf('Error: %s', err.message);
        end

        % Display results
        set(resultText, 'String', resultString);
    end
end

function [theta, time] = projectile_newton_raphson(v0, R, H, g)
    tol = 1e-6;
    max_iter = 100;

    theta = pi / 4;

    for iter = 1:max_iter
        % Function value
        f_theta = (v0 * cos(theta)) * ((v0 * sin(theta) + sqrt((v0 * sin(theta))^2 + 2 * g * H)) / g) - R;

        term1 = -v0 * sin(theta) * ((v0 * sin(theta) + sqrt((v0 * sin(theta))^2 + 2 * g * H)) / g);
        term2 = (v0 * cos(theta)) * (v0 * cos(theta) + g * H / sqrt((v0 * sin(theta))^2 + 2 * g * H)) / g;
        df_theta = term1 + term2;

        if abs(df_theta) < tol
            error('Derivative is too small. Newton-Raphson may fail to converge.');
        end

        theta_next = theta - f_theta / df_theta;

        if abs(theta_next - theta) < tol
            theta = theta_next;
            break;
        end

        theta = theta_next;
    end

    if iter == max_iter
        error('Newton-Raphson did not converge for angle.');
    end

    time = (v0 * sin(theta) + sqrt((v0 * sin(theta))^2 + 2 * g * H)) / g;
end

function plotTrajectory(v0, theta, time, g, ax)
    t = linspace(0, time, 100); % Generate time values
    x = v0 * cos(theta) * t; % Compute x positions
    y = v0 * sin(theta) * t - 0.5 * g * t.^2; % Compute y positions

    axes(ax); % Set the current axes
    cla(ax); % Clear previous plots
    plot(x, y, 'b-', 'LineWidth', 2);
    grid on;
    xlabel('Horizontal Distance (m)');
    ylabel('Vertical Height (m)');
    title('Projectile Trajectory');
    xlim([0, max(x)]);
    ylim([0, max(y) * 1.1]);
end
