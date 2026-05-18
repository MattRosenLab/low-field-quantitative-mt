function MTfit_multi(dataset_list, param_defs, lsfcn, savestr)
%% MTfit_multi
% General simultaneous fitting of multiple datasets with arbitrary parameter sharing.
%
% INPUTS:
%   dataset_list : cell array of dataset structs with fields:
%       .name         : dataset label
%       .model_fun    : @(p, indep_vars{:}) → model output
%       .indep_vars   : cell array of independent variable arrays
%       .data         : dependent variable array
%       .param_names  : cell array of parameter names used by model_fun
%       .weights      : (optional) weight vector for residuals
%       .plot_fun     : (optional) function handle for plotting
%   param_defs : struct defining parameter initial guesses and bounds:
%       param_defs.param_name.init, .lb, .ub
%   lsfcn : string identifying the type of semisolid lineshape to use
%           Options: {'lorentzian','gaussian','superlorentzian','kubo-tomita'}
%           NOTE: 'kubo-tomita' requires sigtau specified in param_defs!
%   savestr    : (optional) string for file name to export results as CSV
%
% Example:
%   MTfit_multi(dataset_list, param_defs, true)

if nargin < 3
    warning('Lineshape function not specified! Setting to Gaussian by default...')
    lsfcn='gaussian';
    savestr = '';
elseif nargin < 4
    savestr = '';
end

%% --- Build parameter vectors from param_defs ---
param_names = fieldnames(param_defs);
x0 = cellfun(@(n) param_defs.(n).init, param_names);
lb = cellfun(@(n) param_defs.(n).lb, param_names);
ub = cellfun(@(n) param_defs.(n).ub, param_names);

%% --- Optimization setup ---
opts = optimoptions('lsqnonlin', ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 16000, ...
    'MaxIterations', 1600, ...
    'UseParallel', false, ...
    'FunctionTolerance', 1e-16, ...
    'StepTolerance', 1e-16, ...
    'OptimalityTolerance', 1e-16);

%% --- Run the fitting ---
objfun = @(x) total_residual(x, param_names, dataset_list, lsfcn);
[x_fit, resnorm, residual, exitflag, output, lambda, jacobian] = ...
    lsqnonlin(objfun, x0, lb, ub, opts);

%% --- Compute confidence intervals ---
ci = nlparci(x_fit, residual, 'jacobian', jacobian);

%% --- Compute biexponential parameters ---
% First, check whether a selective inversion-recovery dataset was specified
name=cell(numel(dataset_list),1);
iv=name;
for ii=1:numel(dataset_list)
    name{ii}=dataset_list{ii}.name;
    iv{ii}=dataset_list{ii}.indep_vars{2};
end

selIRidx=strcmp(name,'Selective IR');
if sum(selIRidx)>0
    selIRflg=true;

    ra       = x_fit(strcmp(param_names,'ra'));
    rb       = x_fit(strcmp(param_names,'rb'));
    mb0      = x_fit(strcmp(param_names,'mb0'));
    r        = x_fit(strcmp(param_names,'r'));
    t2b      = x_fit(strcmp(param_names,'t2b'));
    ma0      = x_fit(strcmp(param_names,'ma0'));
    
    R1p=0.5*(ra + rb + r*mb0 + r + sqrt((ra - rb + r*mb0 - r)^2 ...
        + 4*r*mb0*r));
    R1n=0.5*(ra + rb + r*mb0 + r - sqrt((ra - rb + r*mb0 - r)^2 ...
        + 4*r*mb0*r));
    
    Tp=iv{selIRidx};
    w1p=1/2/Tp*2*pi;

    t2b=t2b*1e-6;   %If t2b is too small (i.e. on order of 1e-5), fitting Jacobian may be singular!

    switch lsfcn %NOTE: there are probably better expressions for 
        % superlorentzian and kubo-tomita
        case {'lorentzian', 'superlorentzian'}
            adjExpFac=1;            %for Lorentzian semisolid pool
        case {'gaussian', 'kubo-tomita'}
            adjExpFac=sqrt(pi/2);   %for Gaussian semisolid pool
    end
    mbstart=exp(-w1p^2 * Tp * t2b * adjExpFac);
    
    Bp=((ma0 - 1)*(ra - R1n) + (ma0 - mbstart)*r*mb0)/(R1p - R1n);
    Bn=-((ma0 - 1)*(ra - R1p) + (ma0 - mbstart)*r*mb0)/(R1p - R1n);

    param_names_biexp = {'R1p','R1n','Bp','Bn','mbstart'};
    bxfit=[R1p,R1n,Bp,Bn,mbstart];
else
    selIRflg=false;
end

%% --- Display results ---
fprintf('\n=== Fitted Parameters ===\n');
for i = 1:numel(param_names)
    fprintf('%-10s : %.5g   (95%% CI: %.5g to %.5g)\n', ...
        param_names{i}, x_fit(i), ci(i,1), ci(i,2));
end

if selIRflg
    fprintf('\n=== Fitted Biexponential Parameters ===\n');
    for i = 1:numel(param_names_biexp)
        fprintf('%-10s : %.5g   \n', ...
            param_names_biexp{i}, bxfit(i));
    end    
end

%% --- Compute R² for each dataset and generate plots ---
fprintf('\n=== R² per dataset ===\n');
p_fit = make_param_struct(x_fit, param_names);
R2_table = table('Size',[numel(dataset_list),2], ...
                 'VariableTypes',{'string','double'}, ...
                 'VariableNames',{'Dataset','R2'});

for d = 1:numel(dataset_list)
    ds = dataset_list{d};
    model_vals = ds.model_fun(p_fit, ds.indep_vars{:}, lsfcn);
    data = ds.data;

    SSres = sum((data(:) - model_vals(:)).^2);
    SStot = sum((data(:) - mean(data(:))).^2);
    R2 = 1 - SSres/SStot;
    fprintf('%-20s : %.4f\n', ds.name, R2);
    R2_table.Dataset(d) = string(ds.name);
    R2_table.R2(d) = R2;

    % Plot if available
    if isfield(ds, 'plot_fun') && ~isempty(ds.plot_fun)
        figure;
        ds.plot_fun(ds.indep_vars, data, model_vals);
        title(sprintf('%s Fit', ds.name));
    end
end

%% --- Export results to CSV ---
if ~isempty(savestr)
%     timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
    filename = sprintf('%s.csv', savestr);

    param_table = table(param_names, x_fit, ci(:,1), ci(:,2), ...
        'VariableNames', {'Parameter','Value','CI_Lower','CI_Upper'});

    if selIRflg
        biexp_table = table(param_names_biexp', bxfit', NaN(5,1), NaN(5,1), ...
            'VariableNames', {'Parameter','Value','CI_Lower','CI_Upper'});
    else
        biexp_table = table;
    end

    R2_rows = [R2_table.Dataset', ...
               arrayfun(@(r) sprintf('%.5f', r), R2_table.R2, 'UniformOutput', false)'];
    R2_rows=reshape(cellstr(R2_rows),numel(dataset_list),[]);
    R2_paramnames = R2_rows(:,1);
    R2_values = cellfun(@str2double, R2_rows(:,2));
    R2_table_export = table(R2_paramnames, R2_values, NaN(size(R2_values)), NaN(size(R2_values)), ...
        'VariableNames', {'Parameter','Value','CI_Lower','CI_Upper'});

    results_table = [param_table; biexp_table; R2_table_export];
    results_table = rows2vars(results_table);    
    writetable(results_table, filename);
    fprintf('Results exported to %s\n', filename);
end

end


%% ------------------------------------------------------------------------
function res = total_residual(x, param_names, dataset_list, lineshape)
% Compute combined residual for all datasets, including optional weighting
p = make_param_struct(x, param_names);
res = [];

for d = 1:numel(dataset_list)
    ds = dataset_list{d};
    model_vals = ds.model_fun(p, ds.indep_vars{:}, lineshape);
    data = ds.data;

    % Base residual
    r = model_vals(:) - data(:);

    % --- Optional weighting ---
    % 1. Per-point weighting
    if isfield(ds, 'weights') && ~isempty(ds.weights)
        if numel(ds.weights) ~= numel(data)
            error('Weight vector for %s has incorrect length.', ds.name);
        end
        r = r .* ds.weights(:);
    end

    % 2. Per-dataset weighting (relative importance)
    if isfield(ds, 'dataset_weight') && ~isempty(ds.dataset_weight)
        r = r * ds.dataset_weight;
    end

    % 3. Normalize to equalize datasets by default (if not overridden)
    nd = numel(data);
    r = r / sqrt(nd);

    % Combine residuals
    res = [res; r];
end
end


%% ------------------------------------------------------------------------
function p = make_param_struct(x, param_names)
% Convert parameter vector into a named struct
p = cell2struct(num2cell(x)', param_names, 2);
end