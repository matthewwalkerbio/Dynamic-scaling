function [Re, Cd] = Re_Cd_measured(data , constants)


%Generates Re and Cd values for an experiment you just did, given your
% measured U and the S of the model and constants specified in the main
% script


% Calculate Re 
Re = constants.rho_f .* data.U .* data.S .* data.real.L  ./ constants.mu;

% Calculate Cd from force balance equation between drag, weight, buoyancy
% V is measured volume of model based on mass
Cd = 2 .* (constants.rho_m - constants.rho_f) .* constants.g .* data.V ./  constants.rho_f  ./ data.U.^2  ./ (data.S.^2 .* data.real.A);


if length(Re) == 1
    
% disp(['Given measured U = ',num2str(data.U,6),' and printed S = ',num2str(data.S,6),',     tank Re = ',num2str(Re,6),'   and tank Cd = ',num2str(Cd,6)]);

end