This code decompose the Stress intesity factors from strain maps directly without the need for integration
Start with generating strain data using the calibration code and then use the output for the main function
The code is self contained and does not need extra functions

M_J_KIII_2D: Calucate the M integral in horizontal and vertical direction. Under development is the L-integral

There are couple of input desks depending on what you want:
* input_desk_Validation: for valdiating the code. It is connected to the "Calibration_2DKIII" function which creates a synthatic field based on the KI-II-III input.
* input_desk_DIC: for 2D and stereo-DIC data
* input_desk_xEBSD: for HR-EBSD data genereated using xEBSD code. for CrossCourt you will need to do the data prepreation yourself, but look to the other inptu desk to see what the format should be and how it is arranged. There is an example data, "Crack_in_Si_XEBSD" but you will need to add  "GetGrainData" to your directory
* Direction_Sweep: this input desk when you want to find the maximum energy direction as the KI-II-II, J and M are claculated from -90 to 90 degrees.

The most critical thing for the code to work is that your data need to be square (number of data in X need to equal number of data in Y) and the discretisation of the data need to be unifrom, i.e., equal spacing

% for this function to work properly, spacing between potins in x and y
% should be the same and the crack should be at the centre (this can be
% done inside this fucntion also)


% this functions accept data from HR-EBSD and sttero-DIC
% The crack needs to be on the centre exactly for the code to
% work and give good results. This code already include the assumption of
% sigma 33 free == 0
% stress in Pa, E in Pa, distance in m
% use variable 'Stiffness' if you are using anistropic material or E and nu
% for istroropic material
% need some unit clbartion specially for mm

% Another option is to input the strain components as a vector matrix with
% 9 columns the first three columns are the x, y and z coordinate in meters. z
% cooridnate can be a zero column. the 4th to the 9th column are the strain
% components arranged as
% Maps = [X(:) Y(:) Z(:) E11(:) E12(:) E13(:) E21(:) E22(:) E23(:) E31(:) E32(:) E33(:)]; 

% if the map is a 2D strain map then zero all out of the plane components

% the material paramters as MatProp.E for Young's Modulus and MatProp.nu
% for Possions ratio. or as a stifness matrix all in Pa
% xEBSD will assume and solve for  plane strain conditions with sigma33 == 0 + no
% volumetric change.
