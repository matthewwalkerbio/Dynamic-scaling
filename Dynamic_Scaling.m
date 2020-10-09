%% Requirements
% 'Statistics and Machine Learning Toolbox'
% Known to work with MATLAB 2016a
%% Set Directory and Folders for saving
% current directory
projectdir = ('C:\\Users\Matth\Documents\\University\\PhD\Writing\\Similitude paper\\Revised Version\\Matlab_Code_Sharing_version\'); % note that to avoid escaped character '\U' use \\ before (e.g.) user 
% projectdir = ('E:\\Similitude paper\\Matlab Code\\Foram Code 12-03-20\\')
cd(projectdir);
% set this folder
addpath(genpath(projectdir));
% Establish a folder for saving the figures to
figures_folder = strcat(projectdir,'figures\');
%% Load data
% Loading most recent file
folder = [cd,'\']; %gives folder to check
searchstr = '*.mat'; % look for files that match pattern, * is wildcard
temp = dir([folder,searchstr]);
files = {temp.name};
datenums = [temp.datenum];%sorts files by time altered
[sorted_datenums, inds] = sort(datenums);
last_datestr = datestr(datenums(inds(end)));
last_file = files{inds(end)};
load([folder,last_file],'mydata'); % only loads variables specified

% Expected data format is as follows, some of these will be created during
% the script. [genus] and [species] are dynamic variables. 
% For sphere data, for calibrating the tank, I used: genus = 'spherical', species = 'sphere'
% this will use a different script to plot some graphs, which are more
% suited to checking the calibration of the tank etc.

% mydata.[genus].[species].real.L % Length Parallel to flow m
%                              .A % Projected Area m^2
%                              .V % Volume of solid material m^3
    % model data, one row per variable:
% mydata.[genus].[species].each_scale.M % mass of model g
%                                    .S % scale of model
%                                    .speed % sinking velocity m s^-1
%                                    .Total_data % 3D sinking speeds
%                                    .Re % Reynolds number
%                                    .Cd % Drag Coefficient
%                                    .Re_sphere % Re equivilent sphere
%                                    .Cd_sphere % Cd equivilent sphere
    % model scales
% mydata.[genus].[species].all_scales % array of scales
    % model sinking velocities
%mydata.[genus].[species].U % one sinking velocity per model scale (RD: is this written?)
    % predicted operating point
%mydata.[genus].[species].p_point.cubic % using a cubic spline
% mydata.[genus].[species]..p_point.linear % using a linear spline
%% Specify Genus and Species (can be looped if required)
genus = 'spherical';
species = 'sphere';

%% Actual sphere sizes used
sphere_sizes = ['10 mm';'12 mm';'14 mm';'18 mm';'20 mm'];

%% Settings

save_figs = false;      % Do you want to save the volume plots?
save_graphs = false;    % Do you want to save the Cd/Re curve?
plot_reps = true;       % Plot all data points and not just the averages?

S_guess = 15;       % Guess for scale factor S that will match real life Re, Cd
Re_cutoff = 260;    % Upper limit for data input to generating solution
Re_prediction_bounds =   [1 60] ;  % [lower  upper] where lower and upper are absolute limits for real life Re 

% Settings for Volume calculations/Plotting 
fancy_V_extrapolation = true; % Extrapolates the volume (based on mass and density of material) of 3D printed particle. Performance may be reduced if turned off.
spline_options.V_S.n_knots = 3;  %number of spline knots for fitted V vs S curve based on mass measurements.  3 is optimal in most cases
spline_options.V_S.knot_placement = 'fixed'; %no effect for n_knots = 2.  for n_knots >= 3, 'free' for more freedom in spline fit, 'fixed' for less freedom
spline_options.V_S.degree = 'linear';   %
 % what units do you want the volume to be displayed in?:
        % unit_value: if 1 then in volume is in m^-3
        %             if 100 then volume is in cm^-3
        %             if 1000 then volume is in mm^-3
        %             if 1000,000 then volume is in µm^-3    
unit_value = 1000;

% Settings for Re/Cd Curve
spline_options.Cd_Re.n_knots = 3; %number of spline knots for fitted Cd vs Re curve.  should possibly be a good amount smaller than number of data points.  and apparently at least 2.
spline_options.Cd_Re.knot_placement = 'free';  %no effect for n_knots = 2.  for n_knots >= 3, 'free' for more freedom in spline fit, 'fixed' for less freedom
spline_options.Cd_Re.concave_up = 'on';  % enforce that spline is concave up everywhere?  'on'  or 'off'
spline_options.Cd_Re.degree = 'cubic';


% Universal Constants
constants.tank.g = 9.8;  % m / s^2
constants.real.g = constants.tank.g;

% Constants in the Tank
constants.tank.tank_dia = 0.8;  % m
constants.tank.rho_m = 1142.705;% 3D printed material density in kg m^-3 
constants.tank.rho_f = 830;     % Tank Fluid density in kg m^-3 
constants.tank.mu = 22E-3 ;     % Tank Fluid viscosity in Pa S

%real life
constants.real.rho_f = 1024; % Seawater Density
constants.real.mu = 1E-3;    % Seawater Viscosity
constants.real.rho_m = 2830; % Particle density (i.e. CaCO3)
%% Manually Remove Outlier data point(s)

% Format example below, genus and species should be specified and 
% numbers are the "bad" scale factors.
% For multiple bad scale factors for a species, can put a vector e.g.
% 'genus1 species1', 12...;
% 'genus1 species2', 12...;
% 'genus2 species1', [12 13 14];
        bads = {'Dentoglobigerina rohri' , 12 ; 
            };
%% obtain the correct data for the species being used

         disp([genus,' ',species]); %print the genus and species
        
        temp = mydata.(genus).(species);% make a temp data storage to work with
        clear data; % clear old data
        
        data.real = temp.real; % set up temp data with real values
        data.U = mydata.(genus).(species).U; % add sinking speed data
        
        vars = {'S','M'}; % get data values S (scale) and M (mass)
        for i = 1:length(vars)% makes vector
            if isfield(temp,'each_scale')
                data.(vars{i})= [temp.each_scale.(vars{i})];
            else
                data.(vars{i})= [temp.data.(vars{i})];
            end
        end
        
        % calculate the volume based on mass and material density
        data.V = data.M / 1000 / constants.tank.rho_m; % model volume (m^3) based on measured mass and assumed density
        %data.V = data.real.V * data.S.^3;  % uncomment this to use scaled volume based on STL instead of measured volume from model mass
        
        % RD: is this needed?
        L_factor = 1;  %what if we defined L differently?  this can be thought of as scaling L by some factor, which will propagate into our Re values
        data.real.L = data.real.L * L_factor;
        
        A_factor = 1;  %what if we defined A differently?  this can be thought of as scaling A by some factor, which will propagate into our Cd values
        data.real.A = data.real.A * A_factor;
        
        % in theory, we can change L_factor and A_factor (e.g. how we define L and A) independently, though we
        % do expect that sensible choices should be dependent in some way
        
        %Check that data.S and Data.U are the same length, if not then skip
        %to end of script
        
        if ~(length(data.U)==length(data.S));
            return
        end
        
        % remove outliery points based on bads list
        if ~isempty(bads)
            [isbad,ind] = ismember([genus,' ',species],bads(:,1));
            if isbad
                fields = {'U','M','V','S'};  %make sure S is the last variable here, otherwise this will fail badly
                for fb = 1:length(fields)
                    data.(fields{fb})( ismember(data.S ,  bads{ind,2}) ) = [];
                end
            end
        end
        
        %RD: delete this?
        % fake data for debugging
        % data.real.L = 0.0003;  data.real.A = 7E-7;  data.real.V = 0.08E-11;
        %         data.S = [10 12 14 17.5 22 35 54];
        %         data.V = data.real.V * data.S.^3;
        %         data.M = data.V * constants.rho_m;
        %         L = data.real.L * data.S;
        %         data.U = 2/9*(constants.rho_m - constants.rho_f) / constants.mu * constants.g * (L/2).^2 ;
        %         data.U = data.S * 1/100 + 0.1;
        
        if spline_options.V_S.n_knots < 3 
            spline_options.V_S.knot_placement = 'fixed'; %for 2 knots
        end
        options = slmset('increasing','on','extrapolation','linear','knots',spline_options.V_S.n_knots,'interiorknots',spline_options.V_S.knot_placement);
        spline_fit_V = slmengine(data.S,data.V,options);
        
        if fancy_V_extrapolation
            clear V_vs_S % clear persistent variables inside function
            predict_V = @(S) V_vs_S(S,spline_fit_V,data);
        else
            predict_V = @(S) slmeval(S, spline_fit_V);
        end
        
        extrapolation_fraction = 0.2;  %just for this graph, code will extrapolate beyond this if needed
        S_plt = linspace(min(data.S)-range(data.S)*extrapolation_fraction,max(data.S)+range(data.S)*extrapolation_fraction,100);
        F_plt = predict_V(S_plt);
  %% Plot measured and expected volume
        [~,inds] = sort(data.S);
        figure(1);  plot(data.S,(data.V*unit_value^3),'bo','markerfacecolor','b');  hold on;
        plot(S_plt,F_plt*unit_value^3,'b--','linewidth',1.5);
            set(gca, 'FontSize', 16)% add for prettiness
        plot(S_plt,(data.real.V .* S_plt.^3 *unit_value^3)  ,'r-','linewidth',1.5);
        grid on;  hold off;
        xlab = xlabel('Scale factor ({\it  S})', 'interpreter', 'tex');
        % check what units the y axis is in and label suitably
        if unit_value == 1;
           ylab = ylabel('Volume ({$\it V$}) (m^3)', 'interpreter', 'tex');
        elseif unit_value == 100;
           ylab = ylabel('Volume ({\it V}) (cm^3)', 'interpreter', 'tex');
        elseif unit_value == 1000;
           ylab = ylabel('Volume ({\it V}) (mm^3)', 'interpreter', 'tex');
        elseif unit_value == 1000000;
           ylab = ylabel('Volume ({\it V}) (µm^3)', 'interpreter', 'tex');
        end 
        set(ylab, 'FontName', 'Cambria Math')
        set(xlab, 'FontName', 'Cambria Math')
        % Specify the legend properties:
        leg = legend({'\it {M}/\it {\rho_{particle}}','{\it H(S)} (Spline Fit)','Scaled up from µCT scan'},'location','northwest', 'interpreter', 'tex'); 
            set(leg, 'FontSize', 12)             % legend font size
            set(leg, 'FontName', 'Cambria Math') % legend font 
        tit_head = title(strcat(genus,{' '},species));  % make the title
            set(tit_head, 'FontAngle', 'italic');       % make the title italic
        if save_figs == 1
            print('-dpng','-r300',[figures_folder,species,' model volume.png']);
        end
 %%           
 % spline_fit_V = [];  % do this if you want to use STL-based
        %predicted V vs S instead of empirical V vs S spline from measured
        %mass
        
        
        % data.U and data.S were loaded from mat file containing Data struct above
        
        [data.Re, data.Cd_tank] = Re_Cd_measured(data, constants.tank);
        
        
        filter = data.Re < Re_cutoff;
        fields = setdiff(fieldnames(data),'real'); %all fields that have values for each expr
        for f = 1:length(fields)
            data.(fields{f}) = data.(fields{f})(filter);
        end
        
        
        % correct tank Cd data to probable Cd in real life using wall effects
        % correction
        data.lambda = data.S * data.real.L  / constants.tank.tank_dia;
        data.Cd_Inf = Cd_tank2Inf(data.Cd_tank, data.Re, data.lambda);
        S = 1;  %real life scale factor always 1
        
        
 %% Cd_Re plot
%  S = 1; % for checking plotting
 Re_cd_fig = figure(2);  % hopefully will be plotted inside predict_real_life.m with data and spline fit 
   if strcmp(genus, 'spherical')  %use special script if plotting sphere data   
        [Re, Cd_Inf, U, realU] = predict_real_life_sphere(data, S, constants.real, Re_prediction_bounds, spline_options.Cd_Re, Re_cd_fig);  %predict Re and Cd for actual foram in either real life or tank, as well as U for real life
                %place the predicted sinking speed, Re and Cd into the dataframe
                %inside predict_real_life in plot_Cd_vs_Re limits on the Re
                %values are imposed, chnage this if the spline is too
                %long/short
   else %other wise just use the standard script
       [Re, Cd_Inf, U, realU] = predict_real_life(data, S, constants.real, Re_prediction_bounds, spline_options.Cd_Re, Re_cd_fig);  %predict Re and Cd for actual foram in either real life or tank, as well as U for real life
   end     
%         disp(['Given input data ',',   actual operating Re = ',num2str(Re,5),'    actual operating  Cd_Inf = ',num2str(Cd_Inf,5),'    Real life U = ',num2str(U,5),'   m/s']);

        % Adding model S prediction
    
    if strcmp(genus, 'spherical'); 
        % special case for the sphere graph
                figure(Re_cd_fig)
         % generate the labels for the graph, manually specified for the
         % sphere dataset
        a = data.S'; b = num2str(a); c = sphere_sizes;%prepares Scale data for the graph labels
        dx = diff(xlim)*1/100; dy = diff(ylim)*1/100; % displacement so the text does not overlay the data points
        text(data.Re+dx, data.Cd_Inf+dy, c,'fontsize',12,'fontweight','bold');%adds labels to the data points
   % goes with figure from predict real life results in above cell
%         title({species,['Predicted operating point:  ','Re = ',num2str(Re),'  Cd_{Inf} = ',num2str(Cd_Inf),'  S = ',num2str(S), ' U = ', num2str(realU)]});
        
% RD: can you check what of this is needed. i know it's not all of it, but
% maybe some (e.g. enforce S in numerical order)
        if plot_reps
            % you sometimes seem to not have an each_scale field....
            if isfield(mydata.(genus).(species) , 'each_scale')
                reps = vertcat( mydata.(genus).(species).each_scale.Total_data ); % row is scale, col is rep
                scales = [mydata.(genus).(species).each_scale.S];  % in case scales are in some weird order
            else
                reps = vertcat( mydata.(genus).(species).data.Total_data ); % row is scale, col is rep
                scales = [mydata.(genus).(species).data.S];  % in case scales are in some weird order
            end
            
            for sc = 1:length(scales)
                for rp = 1:size(reps,2)
                    data_temp = []; data_temp.real = data.real;
                    data_temp.S = scales(sc); data_temp.V = data.V(data.S == scales(sc));
                    data_temp.U = reps(sc,rp);
                    if isempty(data_temp.V) % happens when we've removed this scale via the bads list
                        continue
                    end
                    [Re_temp, Cd_temp] = Re_Cd_measured(data_temp, constants.tank);

                   hold on
                   % plot settings
                   plot(Re_temp,Cd_temp, 'bo','markerfacecolor','cyan','markersize',5)
                   % legend settings
                   lh = legend({'Experimental Data Average','{\it C^E_D(Re)} (Spline Fit)', '{\it C^M_D(Re)} (Equation 11)','Experimental Data Point'},'location','best', 'FontName', 'Cambria Math', 'FontSize', 15, 'interpreter', 'tex');
                   % reorder the legend items for cosmetic reasons
                   neworder = [4, 1, 2, 3]; % the order you want the legend entries to be in:
                                             % data points, data average,
                                             % spline fit, sphere curve
                    lh.PlotChildren = lh.PlotChildren(neworder); % re-order the legend
                    % adjust the x lims if required
%                     xlim([12.5,57.5]);
                   hold off
                    %text(Re_temp,Cd_temp,[num2str(scales(sc)),'-',num2str(rp)],'fontsize',6,'fontweight','bold');   % for figure
                end
            end

            %save the graph
                if save_graphs ==true
                    % png
                    print('-dpng','-r300',[figures_folder,species,' Cd vs Re_',datestr(now, 'dd-mmm-yyyy HH_MM_SS'),'.png']);
                    % matlab figure
                    savefig(Re_cd_fig,[figures_folder,species,' Cd vs Re_',datestr(now, 'dd-mmm-yyyy HH_MM_SS'),'.fig']);
                end   
        end
         
%%    
         %'normal' Re_Cd figure
    else % if not the sphere data then print a normal graph 
        S = 1; % while checking graphing
        figure(Re_cd_fig)
        % also plot new best guess for real life operating point
        hold on
        plot(Re,Cd_Inf,'o','markerfacecolor','r','markeredgecolor','r','markersize',8);
        hold off
        temp = ylim;  ylim([temp(1) max(temp(2),Cd_Inf*1.1)]);  % make sure plot limits include predicted operating point
        temp = xlim;  xlim([min(temp(1),Cd_Inf*0.1)  temp(2)]);  % make sure plot limits include predicted operating point
        %ylim([0,3]); %fix the axes for comparison
        %xlim([1,80]);%
        
%         legend({'Experimental data average','Spline','C_D(Re) from force balance' , 'Ideal sphere' , 'Predicted operating point', 'Experimental data point'},'location','best');
        
        a = data.S'; b = num2str(a); c = cellstr(b); d = strcat('\it S \rm= ', c);%prepares Scale data for the graph labels
        dx = diff(xlim)*1/100; dy = diff(ylim)*1/100; % displacement so the text does not overlay the data points
        text(data.Re+dx, data.Cd_Inf+dy, d,'fontsize',12,'fontweight','bold', 'fontname','Cambria Math', 'interpreter', 'tex');%adds labels to the data points

        % do title below, once we have a predicted S for possible next experiment
        
        % predict nth experiment S and U again given best guess Re and Cd_Inf for real life from predict_real_life above
        
        
        % Re and Cd_Inf should have just come out of predict_real_life.
        % Re is same between real life and tank but Cd_Inf is real life and Cd_tank
        % is tank with wall effects
        do_correction = true;  %do wall effects correction for Cd_tank?  (otherwise Cd_tank = Cd_Inf) Default = true
        
        [S, lambda, Cd_tank, U] = predict_experiment_with_walls(Cd_Inf, Re, data, constants.tank, predict_V, do_correction, S_guess);  %all output is predicted for tank in next experiment
        
        % Title for figure from predict real life results in above cell
        name_species = strcat([genus, ' ', species]);
        title({'\it ',name_species,... % Main title, then subtitle
            ['\rm','Predicted operating point:  ',...
            '{\it Re^O} = ',num2str(round(Re,2)),' ,',... 
              '  {\it C_D^O} = ',num2str(Cd_Inf, '%.2f'),','... %num2str(round(Cd_Inf,2)),' ,',...
              '  {\it S} = ',num2str(round(S,2)),' ,',...
              ' {\it U^O} = ', num2str(round(realU,4)),' ms^{-1}'],... % end subtitle
              },... % end title grouping
          'FontName','Cambria Math', 'interpreter', 'tex');% end title
        % add the individual data points
        if plot_reps
            % you sometimes seem to not have an each_scale field....
            if isfield(mydata.(genus).(species) , 'each_scale')
                reps = vertcat( mydata.(genus).(species).each_scale.Total_data ); % row is scale, col is rep
                scales = [mydata.(genus).(species).each_scale.S];  % in case scales are in some weird order
            else
                reps = vertcat( mydata.(genus).(species).data.Total_data ); % row is scale, col is rep
                scales = [mydata.(genus).(species).data.S];  % in case scales are in some weird order
            end
            
            for sc = 1:length(scales)
                for rp = 1:size(reps,2)
                    data_temp = []; data_temp.real = data.real;
                    data_temp.S = scales(sc); data_temp.V = data.V(data.S == scales(sc));
                    data_temp.U = reps(sc,rp);
                    if isempty(data_temp.V) % happens when we've removed this scale via the bads list
                        continue
                    end
                    [Re_temp, Cd_temp] = Re_Cd_measured(data_temp, constants.tank);
                    % text(Re_temp,Cd_temp,[num2str(scales(sc)),'-',num2str(rp)],'fontsize',6,'fontweight','bold');
                    % % for general use
                   hold on
                   plot(Re_temp,Cd_temp, 'bo','markerfacecolor','cyan','markersize',5)
                     lh = legend({'Experimental Data Average',...   %1
                                  '{\it C^E_D(Re)} (Spline Fit)',... %2
                                  '{\it C_D^F(Re)} Constraint (Equation 7)',...%3
                                  '{\it C^M_D(Re)} (Equation 11)',...%4
                                  'Predicted Operating Point',...   %5
                                  'Experimental Data Point'},...    %6
                                  'location','northeast', 'FontName', 'Cambria Math', 'FontSize', 15,...
                                  'interpreter', 'tex');
                   neworder = [6, 1, 2, 3, 5, 4]; % the order you want the legend entries to be in:
                    lh.PlotChildren = lh.PlotChildren(neworder); % re-order the legend
                   hold off
                    %text(Re_temp,Cd_temp,[num2str(scales(sc)),'-',num2str(rp)],'fontsize',6,'fontweight','bold');   % for figure
                end
            end
            
        end
            %tweak the xlim values for better final figure
%             xlim([20 63])
        
        if save_graphs ==true
            print('-dpng','-r300',[figures_folder,species,' Cd vs Re_',datestr(now, 'dd-mmm-yyyy HH_MM_SS'),'.png']);
            savefig(Re_cd_fig,[figures_folder,species,' Cd vs Re_',datestr(now, 'dd-mmm-yyyy HH_MM_SS'),'.fig']);
        end
    end % end of special case of plotting for sphere vs nonsphere
  
    
        % New data added to data structure 
        mydata.(genus).(species).all_scales.Re = data.Re; %writes out Re # for all of the scales
        mydata.(genus).(species).all_scales.Cd_Inf = data.Cd_Inf; %as above, for Cd_inf 
        mydata.(genus).(species).all_scales.Cd_tank = data.Cd_tank; %Cd in tank