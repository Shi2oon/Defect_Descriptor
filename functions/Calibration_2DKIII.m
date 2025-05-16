function [M2,M4,alldata,Prop,M3] = Calibration_2DKIII(KI,KII,KIII)
% for more details see https://ora.ox.ac.uk/objects/uuid:f2ba08f3-4a27-4619-92ed-bcd3834dadf0/files/d765371972, page 261
%% Input
              close all                   
% Domain size (square, crack tip at centre).
Prop.type         = 'E';    
Prop.units.xy     = 'mm';  % meter (m) or milmeter (mm) or micrometer(um);  
Prop.units.St     = 'Pa'; 
Prop.stressstat   = 'plane_stress'; % 'plane_stress' OR 'plane_strain'
sz = 50;

switch Prop.units.xy
    case 'm'
        saf = 1;
    case 'mm'
        saf = 1e3;
    case 'um'
        saf = 1e6;
    case 'nm'
        saf = 1e9;
end

Prop.E = 210e9;             % Young's Modulus
Prop.nu = 0.3;             % Poisson ratio
G = Prop.E/(2*(1 + Prop.nu));  % Shear modulus
    switch Prop.stressstat
        case 'plane_strain'
            kappa = 3 - (4 .* Prop.nu); % [/]
        case 'plane_stress'
            kappa = (3 - Prop.nu)./(1 + Prop.nu); % [/]
    end                                                           % Bulk modulus                                                  
% SIF loading
KI = KI*1e6;                                                                     % Mode I SIF
KII = KII*1e6;                                                                    % Mode II SIF
KIII = KIII*1e6; 
%{
Maps.Stiffness = [1/Maps.E          -Maps.nu/Maps.E     -Maps.nu/Maps.E 0 0 0
                 -Maps.nu/Maps.E        1/Maps.E        -Maps.nu/Maps.E 0 0 0
                 -Maps.nu/Maps.E    -Maps.nu/Maps.E         1/Maps.E    0 0 0
                  0 0 0     2*(1+Maps.nu)/Maps.E                          0 0
                  0 0 0 0       2*(1+Maps.nu)/Maps.E                        0
                  0 0 0 0 0         2*(1+Maps.nu)/Maps.E];
Maps.Stiffness = Maps.Stiffness^-1;

% Maps.SavingD = [pwd];
% Maps.results = [pwd];
%}
%% Anayltical displacement data.
M2.stepsize = 1/sz*2;
lin = M2.stepsize*(ceil(-1/M2.stepsize)+1/2):M2.stepsize:M2.stepsize*(floor(1/M2.stepsize)-1/2);
[M2.X,M2.Y,M2.Z] = meshgrid(lin,lin,0);

[th,r] = cart2pol(M2.X,M2.Y);
DataSize = [numel(lin),numel(lin),1];
M4.X = M2.X;%*saf;
M4.Y = M2.Y;%;
M4.Z = M2.Z;%;
% displacement data
M4.Ux = ( 0.5*KI/G*sqrt(r/(2*pi)).*(+cos(th/2).*(kappa-cos(th)))+...
              0.5*KII/G*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa+2+cos(th))));
M4.Uy = ( 0.5*KI/G*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa-cos(th)))+...
              0.5*KII/G*sqrt(r/(2*pi)).*(-cos(th/2).*(kappa-2+cos(th))));
M4.Uz = ( 2*KIII/G*sqrt(r/(2*pi)).*sin(th/2));
M4.data = [M4.X(:) M4.Y(:) M4.Z(:) M4.Ux(:) M4.Uy(:) M4.Uz(:)]*saf;

M2.X = M2.X*saf;
M2.Y = M2.Y*saf;
M2.Z = M2.Z*saf;
M4.X = M2.X;%*saf;
M4.Y = M2.Y;%;
M4.Z = M2.Z;%;

% Maps.stepsize = Maps.stepsize*saf;
 
M2.xo = [-0.01;-0.99]*saf;        M2.xm = [0.01;-0.99]*saf;
M2.yo = [0.0026;0.0026]*saf;      M2.ym = [0.03;-0.03]*saf;


%% JMAN approach (without FEM) - Standard J-integral.
% calculating the displacement gradient tensors
[M2.du11,M2.du12,M2.du13] = crackgradient(M4.Ux*saf,M2.stepsize*saf);
[M2.du21,M2.du22,M2.du23] = crackgradient(M4.Uy*saf,M2.stepsize*saf);
[M2.du31,M2.du32,M2.du33] = crackgradient(M4.Uz*saf,M2.stepsize*saf);
alldata = [M2.X(:) M2.Y(:) M2.Z(:) M2.du11(:) M2.du12(:) M2.du13(:)...
     M2.du21(:) M2.du22(:) M2.du23(:) M2.du31(:) M2.du32(:) M2.du33(:)]; 
%% defromation gradient
M3= M2;
            M3.du11 = M3.du11+1;
            M3.du22 = M3.du22+1;
            M3.du33 = M3.du33+1;
M3.data = [M3.X(:) M3.Y(:) M3.Z(:) M3.du11(:) M3.du12(:) M3.du13(:)...
     M3.du21(:) M3.du22(:) M3.du23(:) M3.du31(:) M3.du32(:) M3.du33(:)]; 
%%
% eXX = Maps.E11;
% eYY = Maps.E22;
% eZZ = Maps.E33;
% eXY = 1/2*(Maps.E12+Maps.E21);
% eXZ = 1/2*(Maps.E13+Maps.E31);
% eYZ = 1/2*(Maps.E23+Maps.E32);
% Chauchy stress tensor, assuming linear-elastic, isotropic material
% Maps.S11 = Maps.E/(1-Maps.nu^2)*(eXX+Maps.nu*(eYY+eZZ));
% Maps.S22 = Maps.E/(1-Maps.nu^2)*(eYY+Maps.nu*(eXX+eZZ));
% Maps.S33 = Maps.E/(1-Maps.nu^2)*(eZZ+Maps.nu*(eXX+eYY));
% Maps.S12 = 2*G*eXY;
% Maps.S13 = 2*G*eXZ;
% Maps.S23 = 2*G*eYZ;
end

%% Support function
function [cx,cy,cz]=crackgradient(c,dx)
c=squeeze(c);
[row,~]=size(c);
midr=floor(row/2);
ctop=c(1:midr,:);
cbot=c(midr+1:end,:);
[cxtop,cytop]=gradient(ctop,dx);
[cxbot,cybot]=gradient(cbot,dx);
cx=[cxtop;cxbot];
cy=[cytop;cybot];
cz=zeros(size(cx));
end
