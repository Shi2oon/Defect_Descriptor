function [Maps,alldata]=GetGrainData(fname,Named)
% file name is fname
%% load data
if exist('Named','var') == 0
    [~,Named] = fileparts(fname);
end
Answers = input('Do you want to use a specific grain (S) or the whole map (W), TKD or CrossCourt? ','s');
if strcmpi(Answers, 'S')
    set(0,'DefaultLineMarkerSize',10)
    if ~exist([erase(fname,'.mat') '.mat'],'file')
        fname = erase(fname,'.mat');
        fname = [fname '_XEBSD.mat'];
    end
    load(fname, 'Map_RefID','C_voight', 'GND','Map_EBSD_MTEX','MicroscopeData',...
                    'Grain_Map_A0_sample','Data', 'Grain_Map_rotation_sample', ...
                    'Data_InputMap','Grain_Map_strain_sample','GrainData',...
                    'Grain_Map_stress_sample','Grain_Map_PH2','Grain_Map_MAE2');
    %% plot to select the boundray
    close all;  	warning off;        s3=subplot(1,1,1);
    imagesc(Data_InputMap.X_axis,Data_InputMap.Y_axis,Map_RefID); hold on
    axis off;           axis image;         axis xy;        s3.YDir='reverse';
    colormap jet;       s3.XDir='reverse';
    title('Respond in the Command Line')
    if isempty(GrainData.RefPoint)
        [GrainData.RefPoint] = to_label(Data_InputMap, Map_RefID);
    else
        for i=1:length(GrainData.RefPoint.x)
            GrainData.RefPoint.prop.labels{i} = num2str(i);
        end
    end
    scatter(GrainData.RefPoint.x,GrainData.RefPoint.y,'k','filled');
    scatter(GrainData.RefPoint.x,GrainData.RefPoint.y,'w');
    labelpoints(GrainData.RefPoint.x,GrainData.RefPoint.y,...
        GrainData.RefPoint.prop.labels,'FontSize',20);
    hold off;      set(gcf,'position',[30 50 1300 950])
    pause(0.1)
    Spec    = input('Which Grain you want to explore?   '); % select
    SavingD = fullfile(fileparts(fname),[Named '_Grain no ' num2str(Spec)]);  
    mkdir(SavingD);
    
    %% save data
    A = squeeze(Grain_Map_A0_sample(:,:,Spec,:,:));             % DEFORMAION
    W = squeeze(Grain_Map_rotation_sample(:,:,Spec,:,:));       % ROTATION
    E = squeeze(Grain_Map_strain_sample(:,:,Spec,:,:));         % STRAIN
    S = squeeze(Grain_Map_stress_sample(:,:,Spec,:,:));         % STRESS
    % save GNDS PH and MAE
    [rowC,colZ] = find(S(:,:,1,1));  % linear indices for nonzero element
    Maps.GND    = NaN(size(squeeze(S(:,:,1,1))));
    if ~exist('GND','var'); GND.total = ones(size(squeeze(S(:,:,1,1)))); end
    for i=1:length(colZ);   Maps.GND(rowC(i),colZ(i)) = GND.total(rowC(i),colZ(i)); end
    try Maps.PH  = mean(squeeze(Grain_Map_PH2(:,:,Spec,:)),3);
        Maps.MAE = mean(squeeze(Grain_Map_MAE2(:,:,Spec,:)),3); end
    % Rotation
    Maps.W11 = W(:,:,1,1);      Maps.W12 = W(:,:,1,2);      Maps.W13 = W(:,:,1,3);
    Maps.W21 = W(:,:,2,1);      Maps.W22 = W(:,:,2,2);      Maps.W23 = W(:,:,2,3);
    Maps.W31 = W(:,:,3,1);      Maps.W32 = W(:,:,3,2);      Maps.W33 = W(:,:,3,3);
    % Stress
    Maps.S11 = S(:,:,1,1);      Maps.S12 = S(:,:,1,2);      Maps.S13 = S(:,:,1,3);
    Maps.S21 = S(:,:,2,1);      Maps.S22 = S(:,:,2,2);      Maps.S23 = S(:,:,2,3);
    Maps.S31 = S(:,:,3,1);      Maps.S32 = S(:,:,3,2);      Maps.S33 = S(:,:,3,3);
    % Strain
    Maps.E11 = E(:,:,1,1);      Maps.E12 = E(:,:,1,2);      Maps.E13 = E(:,:,1,3);
    Maps.E21 = E(:,:,2,1);      Maps.E22 = E(:,:,2,2);      Maps.E23 = E(:,:,2,3);
    Maps.E31 = E(:,:,3,1);      Maps.E32 = E(:,:,3,2);      Maps.E33 = E(:,:,3,3);
    % Deformation gradient
    Maps.A11 = A(:,:,1,1);      Maps.A12 = A(:,:,1,2);      Maps.A13 = A(:,:,1,3);
    Maps.A21 = A(:,:,2,1);      Maps.A22 = A(:,:,2,2);      Maps.A23 = A(:,:,2,3);
    Maps.A31 = A(:,:,3,1);      Maps.A32 = A(:,:,3,2);      Maps.A33 = A(:,:,3,3);
    
    % stiffness:  crystal orientation is defined as the rotation that transforms crystal
    % coordinates, i.e., a description of a vector or a tensor with respect to the crystal
    % reference frame, into specimen coordinates, i.e., a description of the same object
    % with respect to a specimen fixed reference frame.
    Maps.R = Map_EBSD_MTEX(sub2ind([MicroscopeData.NROWS,MicroscopeData.NCOLS],...
        GrainData.RefPoint.prop.yi(Spec),GrainData.RefPoint.prop.xi(Spec)))...
        .orientations.matrix;
    Maps.Stiffness = Korsunsky_StiffnessRot(Maps.R,C_voight(:,:,Spec));
    % Dim
    Maps.X   = Data.XSample;    Maps.Y   = Data.YSample;
    Maps.stepsize  =(abs(Maps.X(1,1)-Maps.X(1,2)));
    Maps.Wo  = (1/2).*(Maps.S11.*Maps.E11 + Maps.S12.*Maps.E12 + Maps.S13.*Maps.E13 +...
        Maps.S21.*Maps.E21 + Maps.S22.*Maps.E22 + Maps.S23.*Maps.E23 +...
        Maps.S31.*Maps.E31 + Maps.S32.*Maps.E32 + Maps.S33.*Maps.E33);
    % units (defualt xEBSD units)
    Maps.units.xy = 'um';       Maps.units.S  = 'GPa';  Maps.units.St  = 'GPa';
    Maps.units.E  = 'Abs.';     Maps.units.W = 'rad';
    Map_RefID(Map_RefID~=Spec) = 0;   Maps.RefID = Map_RefID;
    Maps.Mat = Map_EBSD_MTEX.mineralList{Map_EBSD_MTEX.indexedPhasesId(...
       Map_EBSD_MTEX.phase(sub2ind([MicroscopeData.NROWS,MicroscopeData.NCOLS],...
        GrainData.RefPoint.prop.yi(Spec),GrainData.RefPoint.prop.xi(Spec))))};
    
    %
elseif strcmpi(Answers, 'W') || strcmpi(Answers, 'D')
    [Maps]  = loadingXEBSD(fname);
    SavingD = fullfile(fileparts(fname),[Named '_Full_map']);  mkdir(SavingD);
elseif strcmpi(Answers, 'TKD')
    [Maps] = loadingTKD(fname);
    SavingD = fullfile(fileparts(fname),[Named '_TKD']);  mkdir(SavingD);
elseif strcmpi(Answers, 'CrossCourt') || strcmpi(Answers, 'X')
    [Maps] = loadingCrossCourt(fname);
    SavingD = fullfile(fileparts(fname),[Named '_CrossCourt']);  mkdir(SavingD);
end

%% Plot selected
% close all;              s3=subplot(1,1,1);
% imagesc(Data.XSample(1,:),Data.YSample(:,1),squeeze(Grain_Map_stress_sample(:,:,Spec,1,1))); %GPa
% axis image;             set(gca,'Ydir','normal');   %axis off;
% s3.XDir='reverse';      colormap jet;               caxis([-1.5 1.5]);
% c = colorbar;           c.Label.String = 'GPa';     %labelling
% s3.YDir='reverse';      set(gcf,'position',[30 50 1300 950])
% title(['\sigma_{11}^{'  num2str(Spec) '}']);      	xlabel('x[\mum]'); ylabel('y[\mum]');
%   CroppedPlot(Maps,'S')
plotStressandRot(Maps)
saveas(gcf,[SavingD '\Full.tif']); saveas(gcf,[SavingD '\Full.fig']); close 

%% crop data
answer = input('Do you want to crop (C) and/or rotate (R) data (C/R/N)? ','s');
if answer == 'C' || answer == 'c'
    [Maps] = Cropping(Maps,SavingD);
    SavingD = [SavingD '\Cropped Data.mat'];
elseif answer == 'R' || answer == 'r'
    [Maps] = Rot2Crop(Maps,SavingD,Answers);
    SavingD = [SavingD '\Crop & Rot Data.mat'];
else
    SavingD = [SavingD '\Full Data.mat'];
end
plotStressandRot(Maps)
saveas(gcf,[erase(SavingD,'.mat') '.tif']);  
saveas(gcf,[erase(SavingD,'.mat') '.fig']); close 
%{
if length(unique(Maps.RefID))>3;   answer = 1;
	while sum(answer) ~= 0
        contourf(Maps.RefID); axis image; colorbar; colormap jet
        answer = input('Which pseudo grain you want to merge [GB Parent]?, [0 0] to exit ');
        Maps.RefID(Maps.RefID==answer(1))=answer(2);
    end
    Maps.RefID(Maps.RefID==0)=NaN;
    UnIQ = unique(Maps.RefID(:));       UnIQ(isnan(UnIQ))=[];
    for iV=1:length(UnIQ)
        Maps.Stif{iV} = SawpStif(fname);
    end
    Maps.E11(isnan(Maps.RefID))=NaN;
end
%}
%% crack coordinates
%{
% try;    Maps = isNaN_Maps(Maps);       end % trim GB
close all;
imagesc(Maps.X(1,:),Maps.Y(:,1),Maps.E12);
c=colorbar; c.Label.String = ['E_{12} [' Maps.units.xy ']'];
%     caxis([-5e-3 5e-3]);
set(gca,'Ydir','normal');	axis image;colormap jet
title('Answer in the command line');
xlabel(['X [' Maps.units.xy ' ]'],'FontSize',20,'FontName','Times New Roman');
ylabel(['Y [' Maps.units.xy ' ]'],'FontSize',20,'FontName','Times New Roman');
set(gcf,'WindowStyle','normal')
set(gcf,'position',[30 50 1300 950]);
title('E_{12} :: Select the Crack, start from crack tip');
[Maps.xo,Maps.yo] = ginput(2);
title('E_{12} :: Select the Crack mask, start from crack tip');
[Maps.xm,Maps.ym] = ginput(2); close
%}
%% Get stiffness tensor
Maps.SavingD = SavingD;
if strcmpi(Answers, 'w')
    Maps.results = fname;
	% stiffness:  crystal orientation is defined as the rotation that transforms crystal
    % coordinates, i.e., a description of a vector or a tensor with respect to the crystal
    % reference frame, into specimen coordinates, i.e., a description of the same object
    % with respect to a specimen fixed reference frame.
    [Maps.Stiffness,Maps.R] = SawpStif(Maps.results);
end
if strcmpi(answer, 'R')
	% the stiffness tensor is defined with respect to the measurement x and y, once x 
	% once x and y are rotated the same rotation angle need to be applied to the
	% the stiffness tensor to maintain the Euler angles x and y reference.
	% this can be see clearly when the measurement x an y were changed for (001)
	% Single Si Crystal as the [-110] and [110] pole kept moving.
	% See variation with direction in fig. 3 of https://doi.org/10.1109/JMEMS.2009.2039697
    Maps.Stiffness  = S2DRot(Maps.Stiffness,Maps.theta); % rotate the stiffness
end
% [Maps.Crystal] = rotateStrains(Maps,Maps.R);
% plotStressandRot(Maps.Crystal)
% saveas(gcf,[erase(SavingD,'.mat') '_Cry.tif']);  
% saveas(gcf,[erase(SavingD,'.mat') '_Cry.fig']); close

%% save an exit
Maps.type       = 'A';
Maps.stressstat = 'plane_stress';
Maps.Operation  = 'xED';
alldata = [Maps.X(:)        Maps.Y(:)       zeros(size(Maps.X(:))) ...
           Maps.A11(:)-1    Maps.A12(:)     Maps.A13(:)  ...
           Maps.A21(:)      Maps.A22(:)-1   Maps.A23(:)  ...
           Maps.A31(:)      Maps.A32(:)     Maps.A33(:)-1];
if strcmpi(Named, '3D')
        xLin        = Maps.X(1,:);
        [~, index1] = min(abs(xLin-Maps.xo(1)));
        [~, index2] = min(abs(xLin-Maps.xo(2)));
        yLin        = Maps.Y(:,1);
        [~, index3] = min(abs(yLin-Maps.yo(1)));
        Maps.E11(index3,min([index1,index2]):max([index1,index2]))=NaN;
        Maps.E11(Maps.E11==0)=NaN;
        Zo=input('\nWhat is the effective depth of information [nm]?   ');
    alldata = [ Maps.X(:)       Maps.Y(:)           ones(size(Maps.Y(:)))*Zo*1e-3...
                Maps.E11(:)     Maps.E22(:)         Maps.E33(:)...
                Maps.E12(:)     Maps.E13(:)         Maps.E23(:)];
alldata(isnan(alldata(:,4)),:)=[];
elseif strcmpi(Named, '2D')
    alldata = [Maps.X(:)	Maps.Y(:)	Maps.E11(:)	Maps.E22(:)	Maps.E12(:)];
end

save(SavingD,'Maps','alldata'); % save
Maps.SavingD = fileparts(Maps.SavingD);
end
