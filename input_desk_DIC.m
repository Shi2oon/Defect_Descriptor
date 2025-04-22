%% input desk for 3D istropic materials
clc;clear;close all
% Make sure the grid spacing is equal in x and y direction, also make sure
% the crack is at the centre of the map
addpath(genpath([pwd '\functions']));
Prop.E = 210e9;  % Young's Modulus, you can also input a stifness matrix
Prop.nu = 0.3;     % Poisson ratio         
Prop.units.St = 'Pa';  
Prop.units.xy = 'mm';% meter (m) or milmeter (mm) or micrometer(um);
Prop.stressstat = 'plane_stress';% 'plane_stress' OR 'plane_strain'
Prop.Operation = 'DIC'; % 'DIC' for raw 2D and 3D DIC data, 
%                       % 'F' for deformation gradient and 'U' for
%                       displacement gradient
DataDirect = fullfile(pwd,'Data','1KI-2KII-3KII_Data.dat'); % file locatio
% Prop.SavingD = fullfile(pwd,'Data'); % for saving figures and data

% end of user input
Data = importdata(DataDirect);
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(Data.data,Prop);

%% for Ramberg-Osgood tortuous crack
clc;clear;close all
% Make sure the grid spacing is equal in x and y direction, also make sure
% the crack is at the centre of the map
addpath(genpath([pwd '\functions']));
Prop.E = 210e9;  % Young's Modulus, you can also input a stifness matrix
Prop.nu = 0.3;     % Poisson ratio     
Prop.yield = 4E9; % Yield stress in Pa
Prop.Yield_offset = 1.24; % Hardening coefficient
Prop.Exponent = 26.67; % Hardening exponent
Prop.units.St = 'Pa';  
Prop.units.xy = 'm';% meter (m) or milmeter (mm) or micrometer(um);
Prop.stressstat = 'plane_stress';% 'plane_stress' OR 'plane_strain'
Prop.Operation = 'DIC'; % 'DIC' for raw 2D and 3D DIC data
DataDirect = fullfile(pwd,'Data','3160Jm.dat'); % file locatio
% Prop.SavingD = fullfile(pwd,'Data');

% end of user input, cut the contour at 22
Data = importdata(DataDirect);
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(Data.data,Prop);

%% %% input desk for 2D istropic materials
clc;clear;close all
% Make sure the grid spacing is equal in x and y direction, also make sure
% the crack is at the centre of the map
addpath(genpath([pwd '\functions']));
Prop.E = 70e9;  % Young's Modulus, you can also input a stifness matrix
Prop.nu = 0.321;     % Poisson ratio         
Prop.units.St = 'Pa';  
Prop.units.xy = 'nm';% meter (m) or milmeter (mm) or micrometer(um);
Prop.stressstat = 'plane_stress';% 'plane_stress' OR 'plane_strain'
Prop.Operation = 'DIC'; % 'DIC' for raw 2D and 3D DIC data
DataDirect = fullfile(pwd,'Data','Tortuous_Crack_Data.dat'); % file locatio
% Prop.SavingD = fullfile(pwd,'Data','1KI-2KII-3KII_Data');

% end of user input
Data = importdata(DataDirect);
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(Data.data,Prop);