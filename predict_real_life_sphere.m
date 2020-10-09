function [Re, Cd_Inf, U, realU] = predict_real_life_sphere(data, S, constants, bounds, fitting, fig_handle)

% This is a specialised version of predict_real_life for dealing with
% spheres of known size during 'tank calibration'

% After at least 3 experiments, there should be two Cd, Re data points
% and this will interpolate (hopefully) or extrapolate (less hopefully)
% what is happening in real life (Re, Cd, U in real life)
% If extrapolation is still occruing, or there is a large gap in the data
% of the Re/Cd curve then take the output Re  and Cd values and predict the
% next experiment that might be needed to improve the curve (using predict_experiment)



%%

[~,spline_fit] = Cd_vs_Re(NaN, data, fitting); 


if ~isempty(fig_handle)
    figure(fig_handle);
    plot_Cd_vs_Re(data, fitting); hold on;
    xlim([10,60])
    Re_temp = linspace(0.1, max(data.Re)*1.05, 1000); % for graphing

    % Sphere equation
    Cd_sphere = Morrison(Re_temp);
    p = plot(Re_temp, Cd_sphere,'--','linewidth',1.5);  set(p,'color','black');
    xlabel('\it{Re}','FontSize', 20, 'FontName','Cambria Math') % 
    ylabel('\it{C_D}','FontSize', 20, 'FontName','Cambria Math' , 'interpreter', 'tex')% 
    set(gca, 'FontSize', 16)
    hold off
    legend({'Experimental Data Average','{\itC^E_D(Re)} (Spline Fit)', '{\itC^M_D(Re)} (Equation 10)' },'location','best');
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


[~, U, Cd_Inf] = anonfun(Re);
realU = U;
disp(['real Re = ',num2str(Re),'     real Cd_Inf = ',num2str(Cd_Inf),'     real U = ',num2str(U),' m/s'])


