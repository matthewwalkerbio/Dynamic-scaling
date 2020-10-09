function [obj, U, Cd_1, Cd_2] = objective_fun( Re, S, constants, data, spline_fit)
% The output obj should be zero, this is achieved by varying Re and plugging into:
% 1) the definition of Re to get U,
% 2) the force balance equation that yields Cd,
% 3) the empirical relationship between Re and Cd (i.e. a Cd(Re) spline or
% Cd(Re) for a sphere)

% when called by predict_real_life, S = 1




% definition of Re
U = Re ./ constants.rho_f .* constants.mu ./  data.real.L ./   S;

% Cd from force balance
Cd_1 = 2 .* (constants.rho_m - constants.rho_f) .* constants.g .* (S.^3 .* data.real.V) ./  constants.rho_f  ./ U.^2  ./ (S.^2 .* data.real.A);

% Cd from spline fit of Cd vs Re data
Cd_2 = Cd_vs_Re(Re, data, NaN, spline_fit);

obj = Cd_1 - Cd_2; % Find where the two curves intersect i.e. the difference is zero

obj(Re <= 0) = NaN;  % Fix/remove possible negative Re points
