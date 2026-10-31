clc;  % Clear the command window
clear;  % Remove all variables from the workspace
close all;  % Close all figure windows
addpath(genpath(fullfile(pwd, 'functions')));

% Generate synthetic data for a 2D crack analysis
% Replace this block with the experimental HR-EBSD dataset and Prop used
% in the manuscript calculation.
[~, ~, alldata] = Calibration_2DKIII(7.35, 2.15, 0.51);

% Define material properties
Prop.E = 210e9;  % Young's modulus in Pascals
Prop.nu = 0.3;   % Poisson's ratio
Prop.units.St = 'Pa';
Prop.units.xy = 'mm';
Prop.stressstat = 'plane_stress';

% For anisotropic calculations, Prop.Stiffness must contain the 6x6
% stiffness matrix ALREADY expressed in the specimen/measurement frame.
% Example:
% Prop.Stiffness = C_specimen;
%
% IMPORTANT:
% During the direction sweep, the original specimen-frame stiffness is
% transformed independently for each theta. Do not rotate the previously
% rotated stiffness, otherwise rotations accumulate incorrectly.

% Define trial VCE angles
crack_angles = -90:1:90;

% Preallocate outputs
J_values = NaN(length(crack_angles), 1);
KI_values = NaN(length(crack_angles), 1);
KII_values = NaN(length(crack_angles), 1);
KIII_values = NaN(length(crack_angles), 1);
J_STD = NaN(length(crack_angles), 1);
KI_STD = NaN(length(crack_angles), 1);
KII_STD = NaN(length(crack_angles), 1);
KIII_STD = NaN(length(crack_angles), 1);
M_STD = NaN(length(crack_angles), 2);
M_values = NaN(length(crack_angles), 2);
Jv_values = NaN(length(crack_angles), 2);
Jv_STD = NaN(length(crack_angles), 2);
direction = NaN(length(crack_angles), 1);
direction_STD = NaN(length(crack_angles), 1);
Jmax_values = NaN(length(crack_angles), 1);
Jmax_STD = NaN(length(crack_angles), 1);

% Preserve the ORIGINAL material description so each angle is transformed
% from the same specimen-frame state.
Prop0 = Prop;

for i = 1:length(crack_angles)
    theta = crack_angles(i);

    % Same 3D in-plane rotation convention used for the deformation/
    % displacement-gradient tensor in the original toolbox.
    R = [cosd(theta) -sind(theta) 0;
         sind(theta)  cosd(theta) 0;
         0            0           1];

    % ---------------------------------------------------------------------
    % 1) Transform the measured second-order tensor field
    % ---------------------------------------------------------------------
    RotatedAlldata = zeros(size(alldata));

    for idx = 1:size(alldata, 1)
        A = [alldata(idx, 4)  alldata(idx, 5)  alldata(idx, 6);
             alldata(idx, 7)  alldata(idx, 8)  alldata(idx, 9);
             alldata(idx,10)  alldata(idx,11)  alldata(idx,12)];

        % Tensor components in the theta-aligned frame
        A_rot = R * A * R';

        % Keep spatial coordinates unchanged for the VCE-direction sweep.
        % This is therefore NOT a passive-coordinate covariance test.
        RotatedAlldata(idx, 1:3) = alldata(idx, 1:3);
        RotatedAlldata(idx, 4)  = A_rot(1,1);
        RotatedAlldata(idx, 5)  = A_rot(1,2);
        RotatedAlldata(idx, 6)  = A_rot(1,3);
        RotatedAlldata(idx, 7)  = A_rot(2,1);
        RotatedAlldata(idx, 8)  = A_rot(2,2);
        RotatedAlldata(idx, 9)  = A_rot(2,3);
        RotatedAlldata(idx,10)  = A_rot(3,1);
        RotatedAlldata(idx,11)  = A_rot(3,2);
        RotatedAlldata(idx,12)  = A_rot(3,3);
    end

    % ---------------------------------------------------------------------
    % 2) Transform the anisotropic stiffness with THE SAME R
    % ---------------------------------------------------------------------
    Prop_theta = Prop0;

    if isstruct(Prop0) && isfield(Prop0, 'Stiffness')
        % Prop0.Stiffness is assumed to already be in the specimen frame.
        % Korsunsky_StiffnessRot applies the fourth-order/Voigt-equivalent
        % transformation consistent with A_rot = R*A*R'.
        Prop_theta.Stiffness = ...
            Korsunsky_StiffnessRot(R, Prop0.Stiffness);

    elseif isnumeric(Prop0) && isequal(size(Prop0), [6 6])
        % Allow a bare 6x6 stiffness matrix if this input mode is used.
        Prop_theta = Korsunsky_StiffnessRot(R, Prop0);
    end

    % For isotropic E-nu input no stiffness rotation is required because
    % the constitutive tensor is rotationally invariant.

    % ---------------------------------------------------------------------
    % 3) Recompute J and related quantities for this trial VCE direction
    % ---------------------------------------------------------------------
    if i == 1
        [~, KI, KII, KIII, J, M] = ...
            M_J_KIII_2D(RotatedAlldata, Prop_theta);
        oh = length(J.Raw);
    else
        [~, KI, KII, KIII, J, M] = ...
            M_J_KIII_2D(RotatedAlldata, Prop_theta, oh);
    end

    J_values(i) = J.true;
    J_STD(i) = J.div;
    KI_values(i) = KI.true;
    KI_STD(i) = KI.div;
    KII_values(i) = KII.true;
    KII_STD(i) = KII.div;
    KIII_values(i) = KIII.true;
    KIII_STD(i) = KIII.div;
    M_values(i,1:2) = M.true;
    M_STD(i,1:2) = M.div;
    Jv_values(i,1:2) = J.vectorial_true;
    Jv_STD(i,1:2) = J.vectorial_div;
    direction(i) = J.direction_true;
    direction_STD(i) = J.direction_div;
    Jmax_values(i) = J.maxJ_true;
    Jmax_STD(i) = J.maxJ_div;
end

T = table(crack_angles(:), J_values, KI_values, KII_values, KIII_values, ...
    Jv_values(:,1), Jv_values(:,2), M_values(:,1), M_values(:,2), ...
    'VariableNames', ...
    {'Theta_deg','J_J_per_m2','KI_MPa_sqrt_m','KII_MPa_sqrt_m', ...
     'KIII_MPa_sqrt_m','J1_J_per_m2','J2_J_per_m2','M1_J_per_m','M2_J_per_m'});

% Optional:
% writetable(T,'Direction_Sweep_Anisotropic.xlsx');

%% Plot J and K versus VCE direction
set(0,'defaultAxesFontSize',18);
set(0,'DefaultLineMarkerSize',9);

close all;
fig = figure;
set(fig, 'defaultAxesColorOrder', [[0 0 0]; [1 0 0]]);
hold on;

yyaxis left;
errorbar(crack_angles, KI_values, KI_STD, '--ok', ...
    'DisplayName', 'K_I', 'MarkerSize', 11, 'LineWidth', 1.5);
errorbar(crack_angles, KII_values, KII_STD, '-.sk', ...
    'DisplayName', 'K_{II}', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 11, 'LineWidth', 1);
errorbar(crack_angles, KIII_values, KIII_STD, '-<k', ...
    'DisplayName', 'K_{III}', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 11, 'LineWidth', 1);
ylabel('K (MPa m^{0.5})');
hold off;

yyaxis right;
hold on;
errorbar(crack_angles, J_values, J_STD, '--or', ...
    'DisplayName', 'J_{1}^{I+II+III}', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 11, 'LineWidth', 1.5);
errorbar(crack_angles, Jv_values(:,1), Jv_STD(:,1), '-.sr', ...
    'DisplayName', 'J_{1}', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 11, 'LineWidth', 1);
errorbar(crack_angles, Jv_values(:,2), Jv_STD(:,2), '-<r', ...
    'DisplayName', 'J_{2}', 'MarkerSize', 11, 'LineWidth', 1);
hold off;

ylabel('J (J/m^{2})');
xlabel('VCE (deg)');
legend('Location','best');
grid on;
set(gcf, 'position', [30 50 1200 750]);
xticks(-90:15:90);
xlim([-90 90]);

%% Plot M and J vector components versus direction
close all;
fig = figure;
set(fig, 'defaultAxesColorOrder', [[0 0 0]; [1 0 0]]);

yyaxis left;
hold on;
errorbar(crack_angles, M_values(:,1), M_STD(:,1), '-ok', ...
    'DisplayName', 'M_1', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 12, 'LineWidth', 1);
errorbar(crack_angles, M_values(:,2), M_STD(:,2), '-sk', ...
    'DisplayName', 'M_2', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 12, 'LineWidth', 1);
ylabel('M (J/m)');
hold off;

yyaxis right;
hold on;
errorbar(crack_angles, Jv_values(:,1), Jv_STD(:,1), '-^r', ...
    'DisplayName', 'J_1', 'MarkerFaceColor', 'r', ...
    'MarkerSize', 12, 'LineWidth', 1);
errorbar(crack_angles, Jv_values(:,2), Jv_STD(:,2), '-dr', ...
    'DisplayName', 'J_2', 'MarkerSize', 12, 'LineWidth', 1.5);
hold off;

ylabel('J (J/m^{2})');
xlabel('VCE (deg)');
legend('Location','best');
grid on;
set(gcf, 'position', [30 50 1244 643]);

%% J1-J2 trajectory
figure;
plot(Jv_values(:,1), Jv_values(:,2), '.k');
axis equal;
box off;
xlabel('J_1 (J/m^{2})');
ylabel('J_2 (J/m^{2})');
set(gcf, 'position', [30 50 1044 843]);
