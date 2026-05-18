% MTfit_multi_prepInputs.m
%
% This script preps the variables dataset_list and param_defs, which are  
% inputs to function MTfit_multi(), from workspace variables.
%
% This works best if you run it section by section, loading in the
% variables from the .mat file between sections.
%
% !! You will likely need to change the names of the variables in 
% .indep_vars and in .data !!
%
% dataset format:
%       .name         : dataset label
%       .model_fun    : @(p, indep_vars{:}) → model output
%       .indep_vars   : cell array of independent variable arrays
%       .data         : dependent variable array
%       .param_names  : cell array of parameter names used by model_fun
%       .weights      : (optional) weight vector for residuals
%       .plot_fun     : (optional) function handle for plotting
%
%% (1) PREP 1-SIDED Z-SPECTROSCOPY DATA
dataset1 = struct;  %overwrite previous struct variable in workspace
%~~~~~~~~~~~~~~~~~~~~~~~~VALUES TO CHANGE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% dataset1.indep_vars = {mtf_presat_cw.satfrq(:,1), mtf_presat_cw.w1(1,:)};      
% dataset1.data = mtf_presat_cw.specn_sig; 
% dataset1.indep_vars = {[offs_Hz;offs_Hz], w1};      
% dataset1.data(1:52,1:5) = ZSpecMTF30pc_1;
% dataset1.data(53:104,1:5) = ZSpecMTF30pc_2;
dataset1.indep_vars = {offs_Hz, w1};      
dataset1.data = Z_newTRES;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dataset1.name = '1-sided Z-spec';
dataset1.model_fun = @(p, delta, w1, lsfcn) ...
        mt_model_w_dipolar( ...
        p.t2a, p.ra, p.rb, p.mb0, p.r, p.t2b, p.sigtau, p.delta_MT, p.td, ...
        delta, w1, lsfcn);
dataset1.plot_fun = @(iv, data, fit) ...
    plot_zspec(iv{1}, iv{2}, data, fit);

% We need to make sure the .indep_vars are specified as m x 1 and 1 x n,
% respectively
if size(dataset1.indep_vars{1},1)<size(dataset1.indep_vars{1},2)
    dataset1.indep_vars{1}=dataset1.indep_vars{1}';
end
if size(dataset1.indep_vars{2},1)>size(dataset1.indep_vars{2},2)
    dataset1.indep_vars{2}=dataset1.indep_vars{2}';
end

disp('One-sided Z-spectroscopy dataset successfully prepared.')
disp('***REMEMBER TO RE-COMBINE DATASETS INTO VARIABLE dataset_list!!***')


%% (2) PREP 2-SIDED Z-SPECTROSCOPY DATA: ALTERNATING [+ -] FREQUENCIES
dataset2 = struct;  %overwrite previous struct variable in workspace
%~~~~~~~~~~~~~~~~~~~~~~~~VALUES TO CHANGE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dataset2.indep_vars = {mtf_pulsed_double.satfrq(:,1), mtf_pulsed_double.w1(1,:)};
dataset2.data = mtf_pulsed_double.specn_sig;
% dataset2.indep_vars = {offs_Hz_altfreq, w1_altfreq};
% dataset2.data = Z_altfreq_MTF100;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dataset2.name = '2-sided Z-spec, Alternating Frequencies';
dataset2.model_fun = @(p, delta, w1, lsfcn) ...
    mt_model(p.t2a, p.ra, p.rb, p.mb0, p.r, p.t2b, p.sigtau, ...
                   p.delta_MT, delta, w1, lsfcn);
dataset2.plot_fun = @(iv, data, fit) ...
    plot_zspec(iv{1}, iv{2}, data, fit);

% We need to make sure the .indep_vars are specified as m x 1 and 1 x n,
% respectively
if size(dataset2.indep_vars{1},1)<size(dataset2.indep_vars{1},2)
    dataset2.indep_vars{1}=dataset2.indep_vars{1}';
end
if size(dataset2.indep_vars{2},1)>size(dataset2.indep_vars{2},2)
    dataset2.indep_vars{2}=dataset2.indep_vars{2}';
end

disp('Two-sided Z-spectroscopy dataset (alternating frequencies) successfully prepared.')
disp('***REMEMBER TO RE-COMBINE DATASETS INTO VARIABLE dataset_list!!***')


%% (3) PREP 2-SIDED Z-SPECTROSCOPY DATA: SIMULTANEOUS FREQUENCIES
dataset3 = struct;  %overwrite previous struct variable in workspace
%~~~~~~~~~~~~~~~~~~~~~~~~VALUES TO CHANGE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% dataset3.indep_vars = {offs_2freq_Hz, w1adj_2freq_10pct};
dataset3.indep_vars = {offs_Hz_2freq, w1_2freq_adj, 0};
dataset3.data = Z_2freq_newTRES;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dataset3.name = '2-sided Z-spec, Simultaneous Frequencies';
% dataset3.model_fun = @(p, delta, w1) ...
%     mt_model_2freq(p.t2a, p.ra, p.rb, p.mb0, p.r, p.t2b, p.sigtau, ...
%                    p.delta_MT, delta, w1);
dataset3.model_fun = @(p, delta, w1, Txoff, lsfcn) ...
    mt_model_2freq_offFromWater(p.t2a, p.ra, p.rb, p.mb0, p.r, p.t2b, p.sigtau, ...
                   p.delta_MT, delta, w1, Txoff, lsfcn);
dataset3.plot_fun = @(iv, data, fit) ...
    plot_zspec(iv{1}, iv{2}, data, fit);

% We need to make sure the .indep_vars are specified as m x 1 and 1 x n,
% respectively
if size(dataset3.indep_vars{1},1)<size(dataset3.indep_vars{1},2)
    dataset3.indep_vars{1}=dataset3.indep_vars{1}';
end
if size(dataset3.indep_vars{2},1)>size(dataset3.indep_vars{2},2)
    dataset3.indep_vars{2}=dataset3.indep_vars{2}';
end

disp('Two-sided Z-spectroscopy dataset (simultaneous frequencies) successfully prepared.')
disp('***REMEMBER TO RE-COMBINE DATASETS INTO VARIABLE dataset_list!!***')


%% (4) PREP SELECTIVE INVERSION-RECOVERY DATA
dataset4 = struct;  %overwrite previous struct variable in workspace
%~~~~~~~~~~~~~~~~~~~~~~~~VALUES TO CHANGE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% dataset4.indep_vars = {mtf_t1.tau, mtf_t1.pars.p1*1e-6};
% dataset4.data = mtf_t1.norm_sig;
dataset4.indep_vars = {TI_newTRES, Tpuls_newTRES};
dataset4.data = SelIR_newTRES;
% % If you want to weight the residuals towards earlier TI values, uncomment
% % below + change adjfac as desired
% % Weighting based on 1/TI 
% adjfac = 0.05;            % Suggested: 0.01 < adjfac < 0.05   
%                           % 0.02 seems good for MTF 20%
% dataset4.weights = adjfac./dataset4.indep_vars{1} + 1;
% % Weighting based on TI < threshold value
% adjfac = 10;
% TI_threshold = 0.1;       % Threshold TI [s]
% dataset4.weights = adjfac * (dataset4.indep_vars{1} < TI_threshold) + 1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
dataset4.name = 'Selective IR';
dataset4.model_fun = @(p, t1, Tp, lsfcn) ...
    biexpMTfitv2(p.ra, p.rb, p.mb0, p.r, p.ma0, p.t2b, Tp, t1, lsfcn);
dataset4.plot_fun = @(iv, data, fit) ...
    plot_IR(iv{1}, data, fit);

disp('Selective IR dataset successfully prepared.')
disp('***REMEMBER TO RE-COMBINE DATASETS INTO VARIABLE dataset_list!!***')


%% COMBINE ALL DATASETS TOGETHER + CHOOSE LINESHAPE
% Modify the line below if you want to change which datasets are
% simultaneously fitted
dataset_list={dataset3, dataset4};
disp('All datasets combined into cell array dataset_list.')

% lineshape='gaussian';
lineshape='kubo-tomita';
disp(['Semisolid lineshape for modeling: ' lineshape]);

%% SETUP PARAMETER STARTING VALUES AND BOUNDS
% .init     = starting point for fitting
% .lb       = lower bound for fitting
% .ub       = upper bound for fitting
% DK NOTE: Start r at a LOW value (e.g. 1) -- otherwise it might not fit
% properly! For some reason, higher r values barely change when fitting
% both 2-freq Z-spec and selIR
% DK NOTE: "sigtau" is used for Kubo-Tomita fitting (keep below 3!)
param_defs = struct( ...
    'ra',       struct('init',.5,    'lb',0,    'ub',Inf), ...
    'rb',       struct('init',1,    'lb',1,    'ub',1), ...
    'mb0',      struct('init',0.15, 'lb',0,    'ub',0.3), ...
    'r',        struct('init',40,   'lb',0,    'ub',1000), ...
    'ma0',      struct('init',-0.95,'lb',-1,   'ub',-0.8), ...
    't2a',      struct('init',.05,  'lb',0,    'ub',2.5), ... 
    't2b',      struct('init',10,   'lb',0,    'ub',50), ...
    'sigtau',   struct('init',1,    'lb',0.01, 'ub',3), ...
    'delta_MT', struct('init',0,    'lb',0,    'ub',0));
%     'td',       struct('init',10e-3, 'lb',0,    'ub',.3), ...    
%     'sigtau',   struct('init',1,    'lb',1,     'ub',1), ...
%     'delta_MT', struct('init',0,    'lb',-5000,    'ub',5000));    

disp('Input struct variable param_defs successfully initialized.')


%% plotting functions used to prepare dataset structs (DO NOT RUN THIS SECTION)
function plot_zspec(delta, w1, data, fitpts)
    % First, reorder delta and w1, along w/ data and fit
    [delta,sortdeltaidx]=sort(delta,'ascend');
    [w1,sortw1idx]=sort(w1,'ascend');
    data=data(sortdeltaidx,sortw1idx);
    fitpts=fitpts(sortdeltaidx,sortw1idx);
    % Then plot
    colors = lines(numel(w1));
    for i = 1:numel(w1)
%         plot(delta./1000, data(:,i), 'o', 'Color', colors(i,:), ...
%             'DisplayName', sprintf('%s_1 = %.1f rad/s', char(hex2dec('03C9')), w1(i)));
        plot(delta./1000, data(:,i), 'o', 'Color', colors(i,:), ...
            'DisplayName', sprintf('B_1 = %.1f µT', w1(i)/2/pi/42.577));
        hold on;
        plot(delta./1000, fitpts(:,i), '-', 'Color', colors(i,:), 'LineWidth', 2, ...
            'HandleVisibility','off');
    end
    xlabel('Offset [kHz]'); ylabel('Z');
    grid on; axis square;
    title('Z-spectrum');
    legend('show','Location','southeast');
    % Plot so as to avoid M0 offset
    delta_noM0=rmoutliers(delta);
    xmin=min(delta_noM0)/1000;
    xmax=max(delta_noM0)/1000;
    axis square;
    axis([xmin, xmax, 0, 1]);
end

function plot_IR(t, data, fitpts)
    plot(t, data, 'bo'); hold on; 
    plot(t, fitpts, 'r-', 'LineWidth',2);
    xlabel('Inversion time [s]'); ylabel('Signal'); 
    legend('Data','Fit'); 
    grid on; axis square;
end