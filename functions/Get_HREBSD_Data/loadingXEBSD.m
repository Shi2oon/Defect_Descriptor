function [Va] = loadingXEBSD(DirxEBSD)
warning off
try
    DirXEBSD = [erase(DirxEBSD, '.mat') '_XEBSD.mat'];
    load(DirXEBSD,'Maps','GND','iPut','Data','C_voight','Map_RefID','GrainData'); 
catch err
    load(DirxEBSD,'Maps','GND','iPut','Data','C_voight','Map_RefID','GrainData');
end

if exist('iPut','var') || exist('C_voight','var')    	
    if  ~exist('GND','var')
        try
            fname = erase(DirxEBSD, '_XEBSD.mat');
            load([erase(fname,'.ctf') '.mat'],'GND');   
            Va.GND = GND.total; 
        catch
            Va.GND = randi(1e15,size(Maps.W11_F1));
            disp('GNDs does not exisits, a radnom map is created');
        end
    else
        Va.GND     = GND.total;   % from xEBSD  
    end
    % rotation
        Va.W11 = zeros(size(Maps.W11_F1));   Va.W12 = Maps.W12_2;	Va.W13 = Maps.W13_F1;           
        Va.W21 = Maps.W21_F1;   Va.W22 = zeros(size(Maps.W11_F1));	Va.W23 = Maps.W23_F1;    
        Va.W31 = Maps.W31_F1;  	Va.W32 = Maps.W32_F1;   Va.W33 = zeros(size(Maps.W11_F1));  
        Va.PH = Maps.PH_2;      Va.MAE = Maps.MAE_2;

    if  ~exist('C_voight','var')  % old version of xEBSD
        % stress
        Va.S11 = Maps.S11_F;	Va.S12 = Maps.S12_F;    Va.S13 = Maps.S13_F; 
       	Va.S21 = Maps.S12_F;	Va.S22 = Maps.S22_F;    Va.S23 = Maps.S23_F;  
        Va.S31 = Maps.S13_F;	Va.S32 = Maps.S23_F;    Va.S33 = Maps.S33_F;  
        %strain
        Va.E11 = Maps.E11_F;	Va.E12 = Maps.E12_F;    Va.E13 = Maps.E13_F; 
       	Va.E21 = Maps.E12_F;	Va.E22 = Maps.E22_F;    Va.E23 = Maps.E23_F;  
        Va.E31 = Maps.E13_F;	Va.E32 = Maps.E23_F;    Va.E33 = Maps.E33_F; 
        % displacement gradient tensor
        Va.A11 = Va.E11+Va.W11; Va.A12 = Va.E12+Va.W12; Va.A13 = Va.E13+Va.W13;
        Va.A21 = Va.E21+Va.W21; Va.A22 = Va.E22+Va.W22; Va.A23 = Va.E23+Va.W23;
        Va.A31 = Va.E31+Va.W31; Va.A32 = Va.E32+Va.W32; Va.A33 = Va.E33+Va.W33;
        %
        Va.Stiffness  = iPut.stiffnessvalues;
        Va.X   = Data.X;        Va.Y   = Data.Y;
        Va.RefID = Maps.RefID;
        Va.Version = 'xEBSD_V2';
    
    elseif exist('C_voight','var') % NEW xEBSD version
        load(DirxEBSD,'Map_stress_sample','Map_strain_sample','Map_A0_sample');
       Va.RefID = Map_RefID;
        % displacement gradient tensor
Va.A11 = Map_A0_sample(:,:,1,1);        Va.A12 = Map_A0_sample(:,:,1,2);        Va.A13 = Map_A0_sample(:,:,1,3);
Va.A21 = Map_A0_sample(:,:,2,1);        Va.A22 = Map_A0_sample(:,:,2,2);        Va.A23 = Map_A0_sample(:,:,2,3);
Va.A31 = Map_A0_sample(:,:,3,1);        Va.A32 = Map_A0_sample(:,:,3,2);        Va.A33 = Map_A0_sample(:,:,3,3);
        % stress
Va.S11 = Map_stress_sample(:,:,1,1);    Va.S12 = Map_stress_sample(:,:,1,2);    Va.S13 = Map_stress_sample(:,:,1,3);    
Va.S21 = Map_stress_sample(:,:,2,1);    Va.S22 = Map_stress_sample(:,:,2,2);    Va.S23 = Map_stress_sample(:,:,2,3);    
Va.S31 = Map_stress_sample(:,:,3,1);    Va.S32 = Map_stress_sample(:,:,3,2);    Va.S33 = Map_stress_sample(:,:,3,3);    

        % strain
Va.E11 = Map_strain_sample(:,:,1,1);    Va.E12 = Map_strain_sample(:,:,1,2);    Va.E13 = Map_strain_sample(:,:,1,3);
Va.E21 = Map_strain_sample(:,:,2,1);    Va.E22 = Map_strain_sample(:,:,2,2);    Va.E23 = Map_strain_sample(:,:,2,3);
Va.E31 = Map_strain_sample(:,:,3,1);    Va.E32 = Map_strain_sample(:,:,3,2);    Va.E33 = Map_strain_sample(:,:,3,3);
        %
        Va.Stiffness = squeeze(C_voight(:,:,1));   
        Va.X   = Data.XSample;          Va.Y   = Data.YSample;
        if Va.Y(1,1)-Va.Y(1,2)~=0 || Va.X(1,1)-Va.X(2,1)~=0
           [Va.X,Va.Y] = meshgrid(unique(Va.X),unique(Va.Y));
           
        end
        Va.Version = 'xEBSD_V3';
    end   
    % stepsize
    uko = unique(Va.X );                Va.stepsize =(abs(uko(1)-uko(2)));
    
elseif exist('iPut','var') ~= 1 && exist('C_voight','var') ~= 1  && exist('Maps','var') ~= 1
    Va = Maps;
    Va.Version = Maps.Version;
else %if using my code
    [filepath,named,~] = fileparts(DirxEBSD);
    load([filepath '\' erase(named, '_XEBSD') '.mat'],'Maps');
    if isfield(Va,'S21')==0 
        Va.S21 = Va.S12;	 Va.S31 = Va.S13;     Va.S32 = Va.S23;
        Va.E21 = Va.E12;	 Va.E31 = Va.E13;     Va.E32 = Va.E23;
    end
    Va = Maps;
    if isfield(Va,'A11')==0 
        Va.A11 = Va.E11+Va.W11; Va.A12 = Va.E12+Va.W12; Va.A13 = Va.E13+Va.W13;
        Va.A21 = Va.E21+Va.W21; Va.A22 = Va.E22+Va.W22; Va.A23 = Va.E23+Va.W23;
        Va.A31 = Va.E31+Va.W31; Va.A32 = Va.E32+Va.W32; Va.A33 = Va.E33+Va.W33;
    end
    Va.Version = 'xEBSD_V2';
end
Va.Wo = (1/2).*(Va.S11.*Va.E11 + Va.S12.*Va.E12 + Va.S13.*Va.E13 +...
                Va.S21.*Va.E21 + Va.S22.*Va.E22 + Va.S23.*Va.E23 +...
                Va.S31.*Va.E31 + Va.S32.*Va.E32 + Va.S33.*Va.E33);
Va.nu  =  Va.Stiffness(1,2)/(Va.Stiffness(1,1)+ Va.Stiffness(1,2));
Va.E   =  Va.Stiffness(1,1)*(1-2*Va.nu)*(1+Va.nu)/(1-Va.nu);
Va.units.xy = 'um';       Va.units.S  = 'GPa';      Va.units.W = 'rad';
Va.units.E  = 'Abs.';     Va.units.St = 'GPa';
Va.GrainData = GrainData;
end

