function [] = plot_Cd_vs_Re(data, fitting)
% Plot the Re and Cd Data


if isempty(data.Re)  %If there is no data, use the sphere curve
    Re_range = [2 50];
else  % Or if there is data:
    Re_range = [min(data.Re)/1.05   max(data.Re)*1.05]; % bound the limits at 1.1 times the upper and lower Re values
end

% generate a vecotr of Re values across the Re range
Re_vec = linspace(Re_range(1),Re_range(2),100);

% generate a vecotr of Cd values across the Cd range
Cd_vec = Cd_vs_Re(Re_vec, data, fitting);

% Plotting 
plot(data.Re,data.Cd_Inf,'o','markerfacecolor','white','markeredgecolor','blue','markersize',12,'LineWidth',1.5)% Plot the data collected in experiments
hold on
plot(Re_vec, Cd_vec,'b-','linewidth',1.5)  % Plot the theory / fitted curve
hold off
grid on
xlabel('Re');  ylabel('Cd_{Inf}'); % Add labels
legend({'experiments','spline'},'location','best'); % Add legend
% X and Y limits
xlim(Re_range);
ylim([ min(data.Cd_Inf)*0.8 max(data.Cd_Inf)*1.2]);
