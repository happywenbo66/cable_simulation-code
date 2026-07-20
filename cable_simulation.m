% Simulation of the closed steel cable system with variable mass loading
% and nonlinear resistance.
% Features inset plots for t >= 5700 s.
% All plots shifted upward to leave space for a caption at the bottom.
% Author: Wenbo Li

clear; close all; clc;

%% 1. System Parameters
rho = 250;          % Mass linear density (kg/m)
g = 9.81;           % Gravitational acceleration (m/s^2)
angle_deg = 30;     % Inclination angle (degrees)
sin_theta = sind(angle_deg); % Sine of the angle

L12_max = 155;      % Max loaded length of segment 1-2 (m)
L41_max = 55;       % Max loaded length of segment 4-1 (m)
v_inject = 0.01;    % Artificial mass injection rate (m/s)
T_threshold = 3874.95; % Static friction threshold (N)

% Time stepping
dt = 0.01;          % Time step (s)
t_max = 20000;      % Maximum simulation time (s)
t_steady_duration = 100; % Duration to confirm steady state (s)

%% 2. Resistance Force Table
% v_data: Velocity (m/s)
% Fres_data: Resistance Force (N)
v_data = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, ...
          7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5];
Fres_data = [3874.95, 3968.03081, 4283.53906, 4860.47383, 5717.30138, ...
             6830.56669, 8225.56912, 9967.21457, 12036.75162, 14630.12005, ...
             17587.74178, 20747.51237, 24255.48443, 28719.49425, 33762.36005, ...
             38826.12847, 44377.46473, 50429.18701, 56998.5416, 64102.77474, ...
             71759.13271, 80613.68036, 88118.0108, 95974.42929, 104196.22, ...
             112786.9751, 121720.01324, 131025.22385, 141806.92877, 151873.14662];

%% 3. Initialization
L12 = 0;    % Loaded length on segment 1-2 (m)
L41 = 0;    % Loaded length on segment 4-1 (m)
v = 0;      % Linear velocity (m/s)
t = 0;      % Current time (s)

% Pre-allocation for performance
max_steps = round(t_max / dt) + 1;
time_hist = zeros(1, max_steps);
L12_hist = zeros(1, max_steps);
v_hist = zeros(1, max_steps);
a_hist = zeros(1, max_steps);
T_hist = zeros(1, max_steps);
Fres_hist = zeros(1, max_steps);
Fresultant_hist = zeros(1, max_steps);

% Steady state control
steady_counter = 0;
steady_state_reached = false;
idx = 1;

%% 4. Main Simulation Loop
while t < t_max && idx <= max_steps
    
    % --- Update Loaded Lengths ---
    if v > 0 && L41 < L41_max
        L41 = L41 + v * dt;
        L41 = min(L41, L41_max);
    end
    
    if L12 < L12_max
        if L41 < L41_max
            L12 = L12 + v_inject * dt;
        else
            L12 = L12 + v * dt;
        end
        L12 = min(L12, L12_max);
    end
    
    % --- Calculate Forces and Acceleration ---
    T_right = rho * g * sin_theta * L12;
    T_left  = rho * g * sin_theta * L41;
    T = T_right - T_left;
    
    M = rho * (L12 + L41);
    
    if T <= T_threshold
        Fres = T;
        Fresultant = 0;
        a = 0;
        v = 0;
    else
        if v < 0; v = 0; end
        Fres = interp1(v_data, Fres_data, v, 'linear', 'extrap');
        Fresultant = T - Fres;
        a = Fresultant / M;
        v = v + a * dt;
        if v < 0; v = 0; end
    end
    
    % --- Record Histories ---
    time_hist(idx) = t;
    L12_hist(idx)  = L12;
    v_hist(idx)    = v;
    a_hist(idx)    = a;
    T_hist(idx)    = T;
    Fres_hist(idx) = Fres;
    Fresultant_hist(idx) = Fresultant;
    
    % --- Steady State Detection ---
    if abs(a) < 1e-4 && v > 0
        steady_counter = steady_counter + dt;
        if steady_counter >= t_steady_duration
            steady_state_reached = true;
            fprintf('Steady state reached at t = %.2f s with v = %.4f m/s.\n', t, v);
            break;
        end
    else
        steady_counter = 0;
    end
    
    t = t + dt;
    idx = idx + 1;
end

% Trim arrays
idx = idx - 1;
time_hist = time_hist(1:idx);
L12_hist  = L12_hist(1:idx);
v_hist    = v_hist(1:idx);
a_hist    = a_hist(1:idx);
T_hist    = T_hist(1:idx);
Fres_hist = Fres_hist(1:idx);
Fresultant_hist = Fresultant_hist(1:idx);

if ~steady_state_reached
    warning('Steady state detection failed. Simulation forced stopped at 20000 s.');
end

%% 5. Plotting Results with Insets (Adjusted Y-offsets to leave bottom space)
figure('Name', 'Steel Cable Simulation Results', 'Color', 'w', 'Position', [100, 100, 1200, 900]);

% Define time mask for insets (t >= 5700s)
inset_idx = time_hist >= 5700;

% ---------------- Plot 1: L vs Time ----------------
% Bottom Y shifted up from 0.74 to 0.77
ax1 = axes('Position', [0.08, 0.77, 0.65, 0.17]);
plot(time_hist, L12_hist, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('L (m)');
title('Loaded Length L vs Time'); grid on;

if max(time_hist) > 5700
    ax1_inset = axes('Position', [0.76, 0.80, 0.20, 0.14]);
    plot(time_hist(inset_idx), L12_hist(inset_idx), 'b-', 'LineWidth', 1.2);
    xlim([5700, max(time_hist)]);
    xlabel('Time (s)'); ylabel('L (m)');
    title('Zoom (t >= 5700 s)', 'Interpreter', 'none', 'FontSize', 9); grid on;
end

% ---------------- Plot 2: Forces vs Time ----------------
% Bottom Y shifted up from 0.51 to 0.54
ax2 = axes('Position', [0.08, 0.54, 0.65, 0.17]);
plot(time_hist, T_hist, 'r-', 'LineWidth', 1.5); hold on;
plot(time_hist, Fres_hist, 'g-', 'LineWidth', 1.5);
plot(time_hist, Fresultant_hist, 'b--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Force (N)');
title('Forces vs Time'); 
legend('T', 'F_{res}', 'F_{resultant}', 'Location', 'best', 'Interpreter', 'tex'); grid on;

if max(time_hist) > 5700
    ax2_inset = axes('Position', [0.76, 0.57, 0.20, 0.14]);
    plot(time_hist(inset_idx), T_hist(inset_idx), 'r-', 'LineWidth', 1.2); hold on;
    plot(time_hist(inset_idx), Fres_hist(inset_idx), 'g-', 'LineWidth', 1.2);
    plot(time_hist(inset_idx), Fresultant_hist(inset_idx), 'b--', 'LineWidth', 1.2);
    xlim([5700, max(time_hist)]);
    xlabel('Time (s)'); ylabel('Force (N)');
    title('Zoom (t >= 5700 s)', 'Interpreter', 'none', 'FontSize', 9); grid on;
end

% ---------------- Plot 3: Velocity vs Time ----------------
% Bottom Y shifted up from 0.28 to 0.31
ax3 = axes('Position', [0.08, 0.31, 0.65, 0.17]);
plot(time_hist, v_hist, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Velocity v (m/s)');
title('Velocity vs Time'); grid on;

if max(time_hist) > 5700
    ax3_inset = axes('Position', [0.76, 0.34, 0.20, 0.14]);
    plot(time_hist(inset_idx), v_hist(inset_idx), 'k-', 'LineWidth', 1.2);
    xlim([5700, max(time_hist)]);
    xlabel('Time (s)'); ylabel('v (m/s)');
    title('Zoom (t >= 5700 s)', 'Interpreter', 'none', 'FontSize', 9); grid on;
end

% ---------------- Plot 4: Acceleration vs Time ----------------
% Bottom Y shifted up from 0.05 to 0.08, leaving bottom 0.08 space for caption
ax4 = axes('Position', [0.08, 0.08, 0.65, 0.17]);
plot(time_hist, a_hist, 'm-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Acceleration a (m/s^2)');
title('Acceleration vs Time'); grid on;

if max(time_hist) > 5700
    ax4_inset = axes('Position', [0.76, 0.11, 0.20, 0.14]);
    plot(time_hist(inset_idx), a_hist(inset_idx), 'm-', 'LineWidth', 1.2);
    xlim([5700, max(time_hist)]);
    xlabel('Time (s)'); ylabel('a (m/s^2)');
    title('Zoom (t >= 5700 s)', 'Interpreter', 'none', 'FontSize', 9); grid on;
end

%% 6. Space for Annotation / Caption
% Uncomment the code below and replace 'Your caption here...' with your own text:
% text(0.5, 0.04, 'Fig. 3. Simulation results of the steel cable dynamics.', ...
%      'HorizontalAlignment', 'center', 'Units', 'normalized', ...
%      'FontSize', 11, 'Interpreter', 'none');

fprintf('Simulation completed. Total time steps: %d\n', idx);