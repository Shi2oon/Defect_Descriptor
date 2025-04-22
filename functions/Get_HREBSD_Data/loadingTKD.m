function [Va] = loadingTKD(DirxEBSD)
load(DirxEBSD,'Maps','GND','iPut','Data');
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');
load(erase(DirxEBSD,'_XEBSD'),'grains','ebsd');

try
    Va.GND = GND.total;
catch
    Va.GND = randi(1e15,size(Maps.W11_F1));
end
% rotation
Va.W11 = Maps.W11_F1;   Va.W12 = Maps.W12_F1;     Va.W13 = Maps.W13_F1;
Va.W21 = Maps.W21_F1;   Va.W22 = Maps.W22_F1;     Va.W23 = Maps.W23_F1;
Va.W31 = Maps.W31_F1;   Va.W32 = Maps.W32_F1;     Va.W33 = Maps.W33_F1;
Va.PH  = Maps.PH_1;     Va.MAE = Maps.MAE_1;

%strain
Va.E11 = Maps.E11_1;	Va.E12 = Maps.E12_1;    Va.E13 = Maps.E13_1;
Va.E21 = Maps.E12_1;	Va.E22 = Maps.E22_1;    Va.E23 = Maps.E23_1;
Va.E31 = Maps.E13_1;	Va.E32 = Maps.E23_1;    Va.E33 = Maps.E33_1;

% % stress
% Va.S11 = Maps.S11_1;	Va.S12 = Maps.S12_1;    Va.S13 = Maps.S13_1;
% Va.S21 = Maps.S12_1;	Va.S22 = Maps.S22_1;    Va.S23 = Maps.S23_1;
% Va.S31 = Maps.S13_1;	Va.S32 = Maps.S23_1;    Va.S33 = Maps.S33_1;

%new strain vector
plot(grains.boundary);
text(grains,num2str(grains.id));
Spec = input('Which Grain? ');
rot_sample = grains(Spec).meanOrientation.matrix;
C_rotated  = StiffnessRot(rot_sample,iPut.stiffnessvalues);
for ii = 1:size(Maps.E11_1,1)
    for ij = 1:size(Maps.E11_1,2)
        strain_voight=[Va.E11(ii,ij);Va.E22(ii,ij);Va.E33(ii,ij);...
            Va.E12(ii,ij)*2;Va.E13(ii,ij)*2;Va.E23(ii,ij)*2];
        %stress
        stress_voight=C_rotated*strain_voight;
        Va.S11(ii,ij) = stress_voight(1);   Va.S12(ii,ij) = stress_voight(4);
        Va.S13(ii,ij) = stress_voight(5);   Va.S21(ii,ij) = stress_voight(4);
        Va.S22(ii,ij) = stress_voight(2);   Va.S23(ii,ij) = stress_voight(6);
        Va.S31(ii,ij) = stress_voight(5);   Va.S32(ii,ij) = stress_voight(6);
        Va.S33(ii,ij) = stress_voight(3);
    end
end

% displacement gradient tensor
Va.A11 = Va.E11+Va.W11; Va.A12 = Va.E12+Va.W12; Va.A13 = Va.E13+Va.W13;
Va.A21 = Va.E21+Va.W21; Va.A22 = Va.E22+Va.W22; Va.A23 = Va.E23+Va.W23;
Va.A31 = Va.E31+Va.W31; Va.A32 = Va.E32+Va.W32; Va.A33 = Va.E33+Va.W33;
%
Va.Stiffness  = iPut.stiffnessvalues;
Va.X   = Data.X;        Va.Y   = Data.Y;
Va.Version = 'xEBSD_V2';

% stepsize
uko = unique(Va.X );                Va.stepsize =(abs(uko(1)-uko(2)));

Va.Wo = (1/2).*(Va.S11.*Va.E11 + Va.S12.*Va.E12 + Va.S13.*Va.E13 +...
    Va.S21.*Va.E21 + Va.S22.*Va.E22 + Va.S23.*Va.E23 +...
    Va.S31.*Va.E31 + Va.S32.*Va.E32 + Va.S33.*Va.E33);
Va.nu  =  Va.Stiffness(1,2)/(Va.Stiffness(1,1)+ Va.Stiffness(1,2));
Va.E   =  Va.Stiffness(1,1)*(1-2*Va.nu)*(1+Va.nu)/(1-Va.nu);
Va.units.xy = 'um';       Va.units.S  = 'GPa';      Va.units.W = 'rad';
Va.units.E  = 'Abs.';     Va.units.St = 'GPa';
Va.RefID = Maps.NumRoi;
end
