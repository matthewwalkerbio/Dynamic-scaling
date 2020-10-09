%%
default_mu = 22E-3 ; % if temp data is missing, use this
constants.tank.rho_f = 830; %oil


%% load key between genus-species and 3 letter abbrev
% key_filename = 'C:\Hull\Matt\current\Species, size and mass data September 2017_w_High_Res.xlsx';
key_filename = 'H:\University\PhD\Experimental Work\Species, size and mass data September 2017_w_High_Res.xlsx';

[num,txt,raw] = xlsread(key_filename,'File Used');

key.genera = txt(2:end,1);  key.species = txt(2:end,2); key.abbrev = txt(2:end,3);


%% load dates, experiments, and temp

% dates_filename = 'C:\Hull\Matt\current\Dates.xlsx';
dates_filename = 'H:\University\PhD\Laptop MATLAB folder\Sinking Experiments\Dates.xlsx';
[num,txt,raw] = xlsread(dates_filename);

exprs.abbrev = txt(2:end,1);  exprs.scale = num(:,1); exprs.date = txt(2:end,3);  exprs.runs = raw(2:end,4:8);
exprs.datenum = datenum(exprs.date,'dd/mm/yyyy');

ind = find(   cellfun( @isempty, txt(:,10) )      ,1,'first') - 1; %row of last date/temp values
temp_data.date = txt(2:ind,10);  temp_data.temp = num(1:ind-1,10);  temp_data.datenum = datenum(temp_data.date,'dd/mm/yyyy');

figure(234)
plot(temp_data.datenum,temp_data.temp,'o');  datetick;  grid on
temp_spline = pchip(temp_data.datenum,temp_data.temp);
dates_interp = linspace(min(temp_data.datenum),max(temp_data.datenum),1000);
hold on
plot(dates_interp,ppval(temp_spline,dates_interp),'k-');

plot(datenum(exprs.date,'dd/mm/yyyy'),ppval(temp_spline,datenum(exprs.date,'dd/mm/yyyy')),'r.','markersize',5);
hold off
xlabel('date');  ylabel('deg C');


%%
was_ran = cellfun(@isnumeric,exprs.runs);  % true for any time value, false for N/A
ViscData = [];

for r = 1:length(exprs.abbrev)
    abbrev = exprs.abbrev{r};
    scale = exprs.scale(r);
    if isempty(ViscData) || ~isfield(ViscData,abbrev)
        ViscData.(abbrev)(1).scale = scale;
        ViscData.(abbrev)(1).date = NaN(5,1);
    else
        [have_scale,ind] = ismember(scale,[ViscData.(abbrev).scale]);
        if ~have_scale
            ViscData.(abbrev)(end+1).scale = scale;
            ViscData.(abbrev)(end).date = NaN(5,1);
        end
    end
    [~,scale_ind] = ismember(scale,[ViscData.(abbrev).scale]);
    
    
    %     if  any(  was_ran(r,:) & ~isnan(ViscData.(abbrev)(scale_ind).date)  )
    %          error(['Two dates for same species/rep found!']);
    %     end
    
    % NOTE:  by default the code will replace any earlier data with later data
    % in cases of replicates being redone at a later date (technically the last
    % row in the xlsx will be used, which presumably matches the last date....)
    
    ViscData.(abbrev)(scale_ind).date(was_ran(r,:),1) = exprs.datenum(r);
    
end

%%
viscosity_data.temp = [25 10 15 17.5 19.9 22.4]; % deg C
viscosity_data.mu = [14.21 24.83 20.96 20.34 18.41 15.63] * 1E-3;  % cP to Pa s

visc_fit = polyfit(viscosity_data.temp,viscosity_data.mu,1); % 1 for linear fit, I don't think more complicated is warranted without more / less noisy data
temp_plt = linspace(min(viscosity_data.temp),max(viscosity_data.temp),200);
figure(929)
plot(viscosity_data.temp,viscosity_data.mu,'o'); hold on
plot(temp_plt,polyval(visc_fit,temp_plt),'k-');  hold off
grid on
xlabel('deg C');  ylabel('viscosity (Pa s)');

%%
% fill in temps using spline interpolant in case temp wasn't taken on that
% date, and fill in viscosities using linear mu vs temp fit

abbrevs = fieldnames(ViscData);
for a = 1:length(abbrevs)
    for s = 1:length(ViscData.(abbrevs{a}))
        ViscData.(abbrevs{a})(s).temp = ppval(temp_spline,ViscData.(abbrevs{a})(s).date);
        ViscData.(abbrevs{a})(s).mu = polyval(visc_fit,ViscData.(abbrevs{a})(s).temp);
    end
end

%% load data, define constants
% loading most recent file
folder = [cd,'\']; %gives folder to check
searchstr = '*.mat'; % look for files that match pattern, * is wildcard
temp = dir([folder,searchstr]);
files = {temp.name};
datenums = [temp.datenum];%sorts files by time altered
[sorted_datenums, inds] = sort(datenums);
last_datestr = datestr(datenums(inds(end)));
last_file = files{inds(end)};
load([folder,last_file],'mydata'); % only loads variables specified

%% add viscosity to each rep in mydata
genera = fieldnames(mydata);
for g = 1:length(genera) %genera
    species = fieldnames(mydata.(genera{g}));
    for s = 1:length(species) %species
        for sc = 1:length( mydata.(genera{g}).(species{s}).each_scale ) %scales
            
            [genus_ind] = find(ismember(key.genera,genera{g}));
            [species_ind] = find(ismember(key.species,species{s}));
            ind = intersect(genus_ind,species_ind);  %there may be multiple species per genus
            abbrev = key.abbrev{ind};
            if ~isfield(ViscData,abbrev)
                disp([abbrev,' scale ',num2str(mydata.(genera{g}).(species{s}).each_scale(sc).S),' missing from ',dates_filename,'   , using default viscosity']);
                for r = 1:length( mydata.(genera{g}).(species{s}).each_scale(sc).speed )
                    mydata.(genera{g}).(species{s}).each_scale(sc).speed(r).viscosity = default_mu;
                end
                continue
            end
            scales = [ViscData.(abbrev).scale];
            scale_ind = find(scales == mydata.(genera{g}).(species{s}).each_scale(sc).S);
            if isempty(scale_ind)
                disp([abbrev,' scale ',num2str(mydata.(genera{g}).(species{s}).each_scale(sc).S),' missing from ',dates_filename,'   , using default viscosity']);
                for r = 1:length( mydata.(genera{g}).(species{s}).each_scale(sc).speed )
                    mydata.(genera{g}).(species{s}).each_scale(sc).speed(r).viscosity = default_mu;
                end
                continue
            end
            for r = 1:length( mydata.(genera{g}).(species{s}).each_scale(sc).speed )
                mydata.(genera{g}).(species{s}).each_scale(sc).speed(r).viscosity = ViscData.(abbrev)(scale_ind).mu(r);
            end
        end
    end
end
%% add Re to each rep in mydata
genera = fieldnames(mydata);
for g = 1:length(genera) %genera
    species = fieldnames(mydata.(genera{g}));
    for s = 1:length(species) %species
        temp = [];
        for sc = 1:length( mydata.(genera{g}).(species{s}).each_scale ) %scales
           
            for r = 1:length( mydata.(genera{g}).(species{s}).each_scale(sc).speed )
                mydata.(genera{g}).(species{s}).each_scale(sc).speed(r).Re_corr = ...
                    mydata.(genera{g}).(species{s}).real.L * mydata.(genera{g}).(species{s}).each_scale(sc).S * ... model length
                    mydata.(genera{g}).(species{s}).each_scale(sc).speed(r).total * ... speed
                    constants.tank.rho_f  / ... density
                    mydata.(genera{g}).(species{s}).each_scale(sc).speed(r).viscosity;
                temp(sc,r) =  mydata.(genera{g}).(species{s}).each_scale(sc).speed(r).Re_corr;
            end
            if ~isempty(  mydata.(genera{g}).(species{s}).each_scale(sc).speed )
                mydata.(genera{g}).(species{s}).each_scale(sc).Total_Re_corr = [ mydata.(genera{g}).(species{s}).each_scale(sc).speed.Re_corr ];
                mydata.(genera{g}).(species{s}).each_scale(sc).Total_visc = [ mydata.(genera{g}).(species{s}).each_scale(sc).speed.viscosity ];
                mydata.(genera{g}).(species{s}).each_scale(sc).mean_Re_corr = nanmean(mydata.(genera{g}).(species{s}).each_scale(sc).Total_Re_corr);
            else
                mydata.(genera{g}).(species{s}).each_scale(sc).Total_Re_corr = [];
                mydata.(genera{g}).(species{s}).each_scale(sc).Total_visc = [];
                mydata.(genera{g}).(species{s}).each_scale(sc).mean_Re_corr = [];
            end
        end
        temp(temp == 0) = NaN;  % just in case there are zero placeholders
        
        mydata.(genera{g}).(species{s}).all_scales.Re_corr = nanmean(temp,2)';
    end
end

%Save mydata
FileName = ['H:\University\PhD\laptop MATLAB folder\Tracking_&_Data\Figures & mat files','MyData_',datestr(now, 'dd-mmm-yyyy HH_MM_SS'),'.mat'];
save (FileName, 'mydata');

