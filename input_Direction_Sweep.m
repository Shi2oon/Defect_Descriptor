clc;  % Clear the command window
clear;  % Remove all variables from the workspace
close all;  % Close all figure windows
addpath([pwd '\functions'])
% Generate synthetic data for a 2D crack analysis
[~, ~, alldata] = Calibration_2DKIII(7.35, 2.15, 0.51);

% Define material properties
Prop.E = 210e9;  % Young's modulus in Pascals (210 GPa for steel)
Prop.nu = 0.3;  % Poisson's ratio
Prop.units.St = 'Pa';  % Units for stress
Prop.units.xy = 'mm';  % Units for coordinates
Prop.stressstat = 'plane_stress';  % Stress state assumption

% Define the rotation angles for the crack orientation
% crack_angles = -90:5:90;  % Angles from -90 to 90 degrees in increments of 5 degrees
crack_angles = -90:1:90;  % Angles from -90 to 90 degrees in increments of 5 degrees

% Preallocate arrays to store results for the J-integral and stress intensity factors
J_values = NaN(length(crack_angles), 1);  % J-integral values
KI_values = NaN(length(crack_angles), 1);  % Mode I SIF values
KII_values = NaN(length(crack_angles), 1);  % Mode II SIF values
KIII_values = NaN(length(crack_angles), 1);  % Mode III SIF values
J_STD = NaN(length(crack_angles), 1);  % Standard deviation for J-integral
KI_STD = NaN(length(crack_angles), 1);  % Standard deviation for K_I
KII_STD = NaN(length(crack_angles), 1);  % Standard deviation for K_II
KIII_STD = NaN(length(crack_angles), 1);  % Standard deviation for K_III
M_STD = NaN(length(crack_angles), 1);  % Standard deviation for M
M_values = NaN(length(crack_angles), 1);  % Standard deviation for M
Jv_values = NaN(length(crack_angles), 1);  % Standard deviation for J
Jv_STD = NaN(length(crack_angles), 1);  % Standard deviation for J

% Loop through each rotation angle
for i = 1:length(crack_angles)
    theta = crack_angles(i);  % Get the current crack angle
    R = [cosd(theta) -sind(theta) 0;  % Create a rotation matrix
         sind(theta) cosd(theta) 0;
         0 0 1];  % 3D rotation matrix (z-axis rotation)

    % Initialize an array to store rotated data
    RotatedAlldata = zeros(size(alldata));
    
    % Loop through each data point to apply the rotation
    for idx = 1:size(alldata, 1)
        % Extract the strain tensor components from the data
        A = [alldata(idx, 4) alldata(idx, 5) alldata(idx, 6);
             alldata(idx, 7) alldata(idx, 8) alldata(idx, 9);
             alldata(idx, 10) alldata(idx, 11) alldata(idx, 12)];

        % Apply the rotation to the strain tensor
        A_rot = R * A * R';  % Rotate the tensor using the rotation matrix

        % Store the rotated strain components back into the rotated data array
        RotatedAlldata(idx, 1:3) = alldata(idx, 1:3);  % Keep coordinates the same
        RotatedAlldata(idx, 4) = A_rot(1, 1);  % Store rotated strain components
        RotatedAlldata(idx, 5) = A_rot(1, 2);
        RotatedAlldata(idx, 6) = A_rot(1, 3);
        RotatedAlldata(idx, 7) = A_rot(2, 1);
        RotatedAlldata(idx, 8) = A_rot(2, 2);
        RotatedAlldata(idx, 9) = A_rot(2, 3);
        RotatedAlldata(idx, 10) = A_rot(3, 1);
        RotatedAlldata(idx, 11) = A_rot(3, 2);
        RotatedAlldata(idx, 12) = A_rot(3, 3);
    end
    
    % Compute the SIFs and J-integral for the rotated data
    if i == 1  % For the first angle, call the function without previous data
        [~,KI,KII,KIII,J,M] = M_J_KIII_2D(RotatedAlldata, Prop);
        oh = length(J.Raw);  % Store the length of raw J data for subsequent calls
    else  % For subsequent angles, pass the previous data length
        [~,KI,KII,KIII,J,M] = M_J_KIII_2D(RotatedAlldata, Prop, oh);
    end
    
    % Store the computed results in the preallocated arrays
    J_values(i) = J.true;           % Store the true J value
    J_STD(i) = J.div;               % Store the standard deviation for J
    KI_values(i) = KI.true;         % Store the true K_I value
    KI_STD(i) = KI.div;             % Store the standard deviation for K_I
    KII_values(i) = KII.true;       % Store the true K_II value
    KII_STD(i) = KII.div;           % Store the standard deviation for K_II
    KIII_values(i) = KIII.true;     % Store the true K_III value
    KIII_STD(i) = KIII.div;         % Store the standard deviation for K_III
    M_values(i,1:2) = M.true;     % Store the true K_III value
    M_STD(i,1:2) = M.div;         % Store the standard deviation for K_III
    Jv_values(i,1:2) = J.vectorial_true;
    Jv_STD(i,1:2) = J.vectorial_div;
    direction(i) = J.direction_true;
    direction_STD(i) = J.direction_div;
    Jmax_values(i) = J.maxJ_true;           % Store the true J value
    Jmax_STD(i) = J.maxJ_div;               % Store the standard deviation for J
end
T = table(crack_angles(:),J_values,KI_values,KII_values,KIII_values,...
    Jv_values(:,1),Jv_values(:,2),M_values(:,1),M_values(:,2),'VariableNames',...
    {'Theta (deg)','J (J/m2)','KI (MPa m-0.5)','KII (MPa m-0.5)',...
    'KIII (MPa m-0.5)','J1 (J/m2)','J2 (J/m2)','M1 (J/m)','M2 (J/m)'});
% writetable(T,'data.xlsx');

%% Plot the results6
set(0,'defaultAxesFontSize',18);       set(0,'DefaultLineMarkerSize',9)
close all; fig = figure;  % Create a new figure
set(fig, 'defaultAxesColorOrder', [[0 0 0]; [1 0 0]]);  % Set default color order for axes
hold on;  % Hold on to the current plot

% Plot K_I, K_II, and K_III on the left y-axis with error bars
yyaxis left;  % Use the left y-axis
errorbar(crack_angles, KI_values, KI_STD, '--ok', 'DisplayName', 'K_I', 'MarkerSize', 11, 'LineWidth', 1.5);
errorbar(crack_angles, KII_values, KII_STD, '-.sk', 'DisplayName', 'K_{II}', 'MarkerFaceColor', 'k', 'MarkerSize', 11, 'LineWidth', 1);
errorbar(crack_angles, KIII_values, KIII_STD, '-<k', 'DisplayName', 'K_{III}', 'MarkerFaceColor', 'k', 'MarkerSize', 11, 'LineWidth', 1);
ylabel('K (MPa m^{0.5})');  % Label for the left y-axis
hold off;  % Release the hold on the current plot

% Plot J-integral on the right y-axis with error bars
yyaxis right; hold on % Use the right y-axis
errorbar(crack_angles, J_values, J_STD, '--or', 'DisplayName', 'J_{1}^{I+II+III}', 'MarkerFaceColor', 'r','MarkerSize', 11, 'LineWidth', 1.5);
errorbar(crack_angles, Jv_values(:,1), Jv_STD(:,1), '-.sr', 'DisplayName', 'J_{1}', 'MarkerFaceColor', 'r', 'MarkerSize', 11, 'LineWidth', 1);
errorbar(crack_angles, Jv_values(:,2), Jv_STD(:,2), '-<r', 'DisplayName', 'J_{2}', 'MarkerFaceColor', 'r', 'MarkerSize', 11, 'LineWidth', 1);
hold off; ylabel('J (J/m^{2})');  % Label for the right y-axis
% line([])
xlabel('VCE\circ');  % Label for the x-axis (crack angle)
legend('Location', 'best');  % Add a legend to the plot
% title('Stress Intensity Factors and J-integral vs Crack Direction Angles');  % Title for the plot
grid on;  % Enable grid on the plot
set(gcf, 'position', [30 50 1200 750]);  % Set the position and size of the figure window
xticks([-90:15:90]);xlim([-90 90])
% Save the figure in both .fig and .tif formats
[~, objh] = legend('NumColumns',2,'Location','northwest','FontSize',18);
% Instead of "h_legend" use "[~, objh]"
objhl = findobj(objh, 'type', 'line');
set(objhl,'Markersize', 15);
saveas(gcf, 'Ks.fig');  saveas(gcf, 'Ks.tif');  close;  % Close the figure window

%%
close all; fig = figure;  % Create a new figure
set(fig, 'defaultAxesColorOrder', [[0 0 0]; [1 0 0]]);  % Set default color order for axes

% Plot K_I, K_II, and K_III on the left y-axis with error bars
yyaxis left;  % Use the left y-axis
hold on;  % Hold on to the current plot
errorbar(crack_angles, M_values(:,1), M_STD(:,1), '-ok', 'DisplayName', 'M_1', 'MarkerFaceColor', 'k', 'MarkerSize', 12, 'LineWidth', 1);
errorbar(crack_angles, M_values(:,2), M_STD(:,2), '-sk', 'DisplayName', 'M_2', 'MarkerFaceColor', 'k', 'MarkerSize', 12, 'LineWidth', 1);
ylabel('M (J/m)');  % Label for the left y-axis
hold off;  % Release the hold on the current plot

% Plot J-integral on the right y-axis with error bars
yyaxis right;  % Use the right y-axis
hold on;  % Hold on to the current plot
errorbar(crack_angles, Jv_values(:,1), Jv_STD(:,1), '-^r', 'DisplayName', 'J_1', 'MarkerFaceColor', 'r', 'MarkerSize', 12, 'LineWidth', 1);
errorbar(crack_angles, Jv_values(:,2), Jv_STD(:,2), '-dr', 'DisplayName', 'J_2', 'MarkerFaceColor', 'r', 'MarkerSize', 12, 'LineWidth', 1.5);
hold off;  % Hold on to the current plot
ylabel('J (J/m^{2})');  % Label for the right y-axis
xlabel('q\circ');  % Label for the x-axis (crack angle)
legend('Location', 'best');  % Add a legend to the plot
title('M and J-integral vs Crack Direction Angles');  % Title for the plot
grid on;  % Enable grid on the plot
set(gcf, 'position', [30 50 1244 643]);  % Set the position and size of the figure window
% xticks([-90:15:90]);xlim([-90 90])
saveas(gcf, 'Ms.fig');  saveas(gcf, 'Ms.tif');  close

%%
plot(Jv_values(:,1),Jv_values(:,2),'.k'); axis equal; box off
xlabel('J_1 (J/m^{2})')
ylabel('J_2 (J/m^{2})')
set(gcf, 'position', [30 50 1044 843]);
saveas(gcf, 'J1_J2.fig');  saveas(gcf, 'J1_J2.tif');  close