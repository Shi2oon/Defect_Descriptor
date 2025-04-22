% Example, synthetic data
clc;clear;close all
addpath([pwd '\functions'])
[Prop,M4,alldata] = Calibration_2DKIII(3,1,2);
Prop.Operation = 'U'; % for defromation gradient
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(alldata,Prop);
Prop.Operation = 'DIC'; %for DIC data
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(M4.data,Prop);% as desigignated maps

%% M integral verfication after rotation (KI-III values will change)
    theta = J.direction_true;  % Get the current crack angle
    R = [cosd(theta) sind(theta) 0;  % Create a rotation matrix
         -sind(theta) cosd(theta) 0;
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
[~,~,~,~,J2,M2] = M_J_KIII_2D(RotatedAlldata,Prop);

