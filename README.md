# Dynamic-scaling

This is the MATLAB code for the paper to appear in the Journal of Experimental Biology: 
Estimation of sinking velocities using free-falling dynamically scaled models: foraminifera as a test case
by Matthew Walker, Stuart Humphries and Rudi Schuech

The main script is Dynamic_Scaling.m which calls the other functions within this folder. The primary purpose of the code is to predict a scale at which to 3D print an enlarged model to replicate the dynamics of a real organism or particle sinking due to gravity. This is done either using the known Cd(Re) relationship for a sphere or by using a spline fit to existing experimental data collected thus far for enlarged model particles. Existing data needs to include measured sinking velocities, as well as basic parameters of the particle (i.e. length parallel to the flow, projected area and volume of material) and the fluid (density, viscosity). All values should be given in SI units.  The code also predicts the sinking velocity of the real particle, given the experimentally observed sinking velocities of a number of scaled models.

In the first instance contact either Matthew Walker (matthewwalkerbio@gmail.com) or Rudi Schuech (rudi.schuech@gmail.com) for assistance. Examplar data avaliable on request. 


The data format that we used is as follows:
% Expected data format is as follows, some of these will be created during
% the script. [genus] and [species] are dynamic variables. 
% For sphere data, for calibrating the tank, I used: genus = 'spherical', species = 'sphere'
% this will use a different script to plot some graphs, which are more
% suited to checking the calibration of the tank etc.

% mydata.[genus].[species].real.L % Length Parallel to flow m
% mydata.[genus].[species].real.A % Projected Area m^2
% mydata.[genus].[species].real.V % Volume of solid material m^3
    % model data, one row per variable:
% mydata.[genus].[species].each_scale.M % mass of model g
% mydata.[genus].[species].each_scale.S % scale of model
% mydata.[genus].[species].each_scale.speed % sinking velocity m s^-1
% mydata.[genus].[species].each_scale.Total_data % 3D sinking speeds
% mydata.[genus].[species].each_scale.Re % Reynolds number
% mydata.[genus].[species].each_scale.Cd % Drag Coefficient
% mydata.[genus].[species].each_scale.Re_sphere % Re equivilent sphere
% mydata.[genus].[species].each_scale.Cd_sphere % Cd equivilent sphere
    % model scales
% mydata.[genus].[species].all_scales % array of scales
    % model sinking velocities
%mydata.[genus].[species].U % one sinking velocity per model scale (RD: is this written?)
    % predicted operating point
%mydata.[genus].[species].p_point.cubic % using a cubic spline
% mydata.[genus].[species]..p_point.linear % using a linear spline
