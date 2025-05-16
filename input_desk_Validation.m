% Example, synthetic displacement gradient list data
clc;clear;close all
addpath([pwd '\functions'])
[~,~,alldata,Prop] = Calibration_2DKIII(3,1,2);
Prop.Operation = 'U'; % for displacement gradient
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(alldata,Prop);
%% displacement gradient organised data into matrix
clc;clear;close all
addpath([pwd '\functions'])
[M2,~,~,Prop] = Calibration_2DKIII(3,1,2);
Prop.Operation = 'U'; % for displacement gradient
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(M2,Prop);
%% displacement organised data into list
clc;clear;close all
addpath([pwd '\functions'])
[~,M4,~,Prop] = Calibration_2DKIII(3,1,2);
Prop.Operation = 'DIC'; %for DIC data
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(M4.data,Prop);
%% displacement organised data into matrix
clc;clear;close all
addpath([pwd '\functions'])
[~,M4,~,Prop] = Calibration_2DKIII(3,1,2);
Prop.Operation = 'DIC'; %for DIC data
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(M4,Prop);

%% Example, synthetic deformation gradient list data
clc;clear;close all
addpath([pwd '\functions'])
[~,~,~,Prop,M3] = Calibration_2DKIII(3,1,2);
Prop.Operation = 'F'; % for deformation gradient
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(M3.data,Prop);
%% deformation gradient organised data into matrix
clc;clear;close all
addpath([pwd '\functions'])
[~,~,~,Prop,M3] = Calibration_2DKIII(3,1,2);
Prop.Operation = 'F'; % for deformation gradient
[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(M3,Prop);

