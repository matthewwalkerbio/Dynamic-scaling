function [S, lambda, Cd_tank, U] = predict_experiment_with_walls(Cd_Inf, Re, data, constants, predict_V, do_correction, S_guess)

% wrapper for orig predict_experiment that takes in Cd_Inf and Re (both our
% best guess for real life) and the constants and do_correction which if true includes wall effects in calculations, if false ignores them, 
% and outputs
% S, scale factor for next experiment
% lambda, wall effects parameter for next experiment
% Cd_tank, predicted Cd in tank in next experiment
% U, predicted speed for next experiment

if isempty(predict_V)
    disp('No empirical spline available for V vs S; falling back to expected V vs S from scaling up STL');
end

obj_fun = @(S_temp) obj_predict_experiment(Cd_Inf, Re, S_temp, data, constants, predict_V,  do_correction);  %keep everything but S_guess constant in search for solution

S = fzero(obj_fun, S_guess);

[~, lambda, Cd_tank, U] = obj_predict_experiment(Cd_Inf, Re, S, data, constants, predict_V, do_correction);  % get out final lambda, Cd_tank, U going with MATLAB's S solution

% disp(['Given input Re = ',num2str(Re,5),' and Cd_Inf = ',num2str(Cd_Inf,5),', suggested scale factor S = ',num2str(S,5),'   predicted tank speed U = ',num2str(U,5),' m/s','   predicted Cd_tank = ',num2str(Cd_tank,5),'   wall effects lambda = ',num2str(lambda,5)]);
