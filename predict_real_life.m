function [Re, Cd_Inf, U, realU] = predict_real_life(data, S, constants, bounds, fitting, fig_handle)


% After at least 3 experiments, there should be two Cd, Re data points
% and this will interpolate (hopefully) or extrapolate (less hopefully)
% what is happening in real life (Re, Cd, U in real life)
% If extrapolation is still occruing, or there is a large gap in the data
% of the Re/Cd curve then take the output Re  and Cd values and predict the
% next experiment that might be needed to improve the curve (using predict_experiment)


%%

[~,spline_fit] = Cd_vs_Re(NaN, data, fitting); % A initial run to save the polyfit coeffs from fitting the data


% plot spline fit right away, before real life prediction calculation
% (which may fail)
if ~isempty(fig_handle)
    figure(fig_handle);
    plot_Cd_vs_Re(data, fitting); hold on; % Re prediction limits in here
    Re_temp = linspace(min(data.Re)*0.5, max(data.Re)*2, 1000); % min and max values for the spline fitting, along with number of sample points

    % Cd(Re) from force balance
    % operating point visually is the intersection between this relationship between
    % Cd(Re) and the actual Cd(Re) spline curve
    % S is always input as 1 here for real life
    Cd_force_balance = 2 .* (constants.rho_m - constants.rho_f) .* constants.g .* (S.^3 .* data.real.V) ./  constants.rho_f  ./ ( Re_temp ./ constants.rho_f .* constants.mu ./  data.real.L ./   S ).^2  ./ (S.^2 .* data.real.A);
    
    plot(Re_temp, Cd_force_balance,'k:','linewidth',1.5);

    % For comparison plot Morrison's (2013) sphere equation
    Cd_sphere = Morrison(Re_temp);
    p = plot(Re_temp, Cd_sphere,'--','linewidth',1.5);  set(p,'color',repmat(0.3,1,3));
    xlabel('\it{Re}','FontSize', 20, 'FontName','Cambria Math') % X Label
    ylabel('\it{C_D}','FontSize', 20, 'FontName','Cambria Math', 'interpreter', 'tex')% Y label
    set(gca, 'FontSize', 16)% overall font size
    hold off
    legend({'Experimental Data Average','{\itC^E_D(Re)} (Spline Fit)', '{\itC^M_D(Re)} (Equation 10)' },'location','best'); % add legend
    drawnow;
end


anonfun = @(x)objective_fun(x(1), S, constants, data, spline_fit);  %in the quest for the best answer, fmincon or patternsearch is allowed to vary  Re but nothing else
num_guesses = 10;
guesses = logspace(log10(bounds(1)),log10(bounds(2)),num_guesses);
options = optimset('display','off');


Re = NaN(1,length(bounds));
for g = 1:length(guesses)
    
    Re(g) = fzero(anonfun,guesses(g),options);
end

Re = min(Re);  % assume correct solution is always the lowest Re intersection unless proven otherwise

    
    
% Re = fzero(anonfun,bounds(2));
% 
% while Re < bounds(1) || Re > bounds(3) % fzero may end up at a bad point sometimes
%     guess = bounds(1) + (bounds(3) - bounds(1))*rand; %random guess between the bounds
%     Re = fzero(anonfun,guess);
% end

[~, U, Cd_Inf] = anonfun(Re);
realU = U;
disp(['real Re = ',num2str(Re),'     real Cd_Inf = ',num2str(Cd_Inf),'     real U = ',num2str(U),' m/s'])


