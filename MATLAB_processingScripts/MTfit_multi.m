function MTfit_multi(dataset_list, param_defs, pars, savestr)
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
%   pars : struct with fitting and CI options:
%       .lsfcn          : string identifying the semisolid lineshape to use
%                         Options: {'lorentzian','gaussian','superlorentzian','kubo-tomita'}
%                         NOTE: 'kubo-tomita' requires sigtau specified in param_defs!
%       .ci_method       : 'nlparci' (default), 'montecarlo', or 'profile'
%       .mc_n_iter       : number of Monte Carlo iterations (default: 500)
%       .mc_use_parallel : use parfor for MC refit loop (default: false)
%       .pl_n_grid       : profile likelihood grid points per side of optimum (default: 50)
%       .pl_range_scale  : grid half-width as multiple of nlparci half-width (default: 4)
%       .pl_alpha        : significance level for profile CI threshold (default: 0.05)
%       .pl_show_profiles: show diagnostic RSS profile plot per parameter (default: false)
%   savestr    : (optional) string for file name to export results as CSV
%
% Example:
%   pars.lsfcn     = 'gaussian';
%   pars.ci_method = 'montecarlo';
%   pars.mc_n_iter = 1000;
%   MTfit_multi(dataset_list, param_defs, pars)

if nargin < 3 || isempty(pars)
    pars = struct();
end
if nargin < 4
    savestr = '';
end

%% --- Parse pars struct with defaults ---
if ~isfield(pars, 'lsfcn') || isempty(pars.lsfcn)
    warning('Lineshape function not specified in pars.lsfcn! Setting to Gaussian by default...')
    pars.lsfcn = 'gaussian';
end
if ~isfield(pars, 'ci_method'),       pars.ci_method       = 'nlparci'; end
if ~isfield(pars, 'mc_n_iter'),       pars.mc_n_iter       = 500;       end
if ~isfield(pars, 'mc_use_parallel'),  pars.mc_use_parallel  = false;     end
if ~isfield(pars, 'pl_n_grid'),        pars.pl_n_grid        = 50;        end
if ~isfield(pars, 'pl_range_scale'),   pars.pl_range_scale   = 4;         end
if ~isfield(pars, 'pl_alpha'),         pars.pl_alpha         = 0.05;      end
if ~isfield(pars, 'pl_show_profiles'), pars.pl_show_profiles = false;     end

lsfcn = pars.lsfcn;

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
use_mc = strcmpi(pars.ci_method, 'montecarlo');
use_pl = strcmpi(pars.ci_method, 'profile');

if use_pl
    % Profile likelihood CIs:
    %   For each parameter k, fix it at a grid of values sweeping outward
    %   from the optimum and re-optimise all other parameters. The CI spans
    %   the range where the profile RSS stays within the F(1, n-p) threshold.
    %   Grid range is initialised from nlparci half-widths × pl_range_scale.
    fprintf('Running profile likelihood CI estimation...\n');

    ci_nlp_internal = nlparci(x_fit, residual, 'jacobian', jacobian);

    n_params  = numel(param_names);
    n_eff     = numel(residual);
    fcrit     = finv(1 - pars.pl_alpha, 1, n_eff - n_params);
    rss_thresh = resnorm * (1 + fcrit / (n_eff - n_params));

    opts_pl = optimoptions('lsqnonlin', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 16000, ...
        'MaxIterations', 1600, ...
        'UseParallel', false, ...
        'FunctionTolerance', 1e-16, ...
        'StepTolerance', 1e-16, ...
        'OptimalityTolerance', 1e-16);

    ci = NaN(n_params, 2);

    for k = 1:n_params
        % Skip fixed parameters (lb == ub — no free variation to profile)
        if lb(k) == ub(k)
            fprintf('  Skipping %d/%d: %s (fixed parameter)\n', k, n_params, param_names{k});
            continue;
        end

        fprintf('  Profiling %d/%d: %s\n', k, n_params, param_names{k});

        % Grid half-width from scaled nlparci estimate
        half_w = pars.pl_range_scale * max(x_fit(k) - ci_nlp_internal(k,1), ...
                                           ci_nlp_internal(k,2) - x_fit(k));
        if ~isfinite(half_w) || half_w == 0
            half_w = 0.5 * max(abs(x_fit(k)), 1e-6);
        end
        lo_k = max(lb(k), x_fit(k) - half_w);
        hi_k = min(ub(k), x_fit(k) + half_w);

        % Reduced bounds (all parameters except k)
        lb_red = lb([1:k-1, k+1:end]);
        ub_red = ub([1:k-1, k+1:end]);

        % Sweep rightward from optimum, warm-starting each step
        right_grid = linspace(x_fit(k), hi_k, pars.pl_n_grid + 1)';
        right_rss  = zeros(pars.pl_n_grid + 1, 1);
        right_rss(1) = resnorm;
        x0_r = x_fit([1:k-1, k+1:end]);
        n_right = 1;
        for g = 2:numel(right_grid)
            v = right_grid(g);
            objfun_pl = @(xr) total_residual([xr(1:k-1); v; xr(k:end)], ...
                param_names, dataset_list, lsfcn);
            x_red_opt = lsqnonlin(objfun_pl, x0_r, lb_red, ub_red, opts_pl);
            r_pl = total_residual([x_red_opt(1:k-1); v; x_red_opt(k:end)], ...
                param_names, dataset_list, lsfcn);
            right_rss(g) = sum(r_pl.^2);
            x0_r = x_red_opt;
            n_right = g;
            if right_rss(g) > rss_thresh; break; end
        end
        right_grid = right_grid(1:n_right);
        right_rss  = right_rss(1:n_right);

        % Sweep leftward from optimum, warm-starting each step
        left_grid = linspace(x_fit(k), lo_k, pars.pl_n_grid + 1)';
        left_rss  = zeros(pars.pl_n_grid + 1, 1);
        left_rss(1) = resnorm;
        x0_l = x_fit([1:k-1, k+1:end]);
        n_left = 1;
        for g = 2:numel(left_grid)
            v = left_grid(g);
            objfun_pl = @(xr) total_residual([xr(1:k-1); v; xr(k:end)], ...
                param_names, dataset_list, lsfcn);
            x_red_opt = lsqnonlin(objfun_pl, x0_l, lb_red, ub_red, opts_pl);
            r_pl = total_residual([x_red_opt(1:k-1); v; x_red_opt(k:end)], ...
                param_names, dataset_list, lsfcn);
            left_rss(g) = sum(r_pl.^2);
            x0_l = x_red_opt;
            n_left = g;
            if left_rss(g) > rss_thresh; break; end
        end
        left_grid = left_grid(1:n_left);
        left_rss  = left_rss(1:n_left);

        % Upper CI: interpolate first threshold crossing going right
        hi_cross = find(right_rss(1:end-1) <= rss_thresh & right_rss(2:end) > rss_thresh, 1);
        if isempty(hi_cross)
            warning('Profile: upper CI for ''%s'' not resolved; reporting grid boundary.', ...
                param_names{k});
            ci(k,2) = hi_k;
        else
            t = (rss_thresh - right_rss(hi_cross)) / ...
                (right_rss(hi_cross+1) - right_rss(hi_cross));
            ci(k,2) = right_grid(hi_cross) + t*(right_grid(hi_cross+1) - right_grid(hi_cross));
        end

        % Lower CI: interpolate first threshold crossing going left
        lo_cross = find(left_rss(1:end-1) <= rss_thresh & left_rss(2:end) > rss_thresh, 1);
        if isempty(lo_cross)
            warning('Profile: lower CI for ''%s'' not resolved; reporting grid boundary.', ...
                param_names{k});
            ci(k,1) = lo_k;
        else
            t = (rss_thresh - left_rss(lo_cross)) / ...
                (left_rss(lo_cross+1) - left_rss(lo_cross));
            ci(k,1) = left_grid(lo_cross) + t*(left_grid(lo_cross+1) - left_grid(lo_cross));
        end

        % Optional diagnostic profile plot
        if pars.pl_show_profiles
            full_grid = [flip(left_grid(2:end)); right_grid];
            full_rss  = [flip(left_rss(2:end));  right_rss];
            figure;
            plot(full_grid, full_rss, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8);
            hold on;
            yline(rss_thresh, 'r--', 'LineWidth', 1.5, 'Label', '95% threshold');
            xline(x_fit(k), 'k:', 'LineWidth', 1.2);
            if isfinite(ci(k,1)), xline(ci(k,1), 'r-', 'LineWidth', 1.2); end
            if isfinite(ci(k,2)), xline(ci(k,2), 'r-', 'LineWidth', 1.2); end
            xlabel(param_names{k}, 'Interpreter', 'none');
            ylabel('Profile RSS');
            title(sprintf('Profile likelihood: %s', param_names{k}), 'Interpreter', 'none');
            hold off;
        end
    end

elseif use_mc
    % Monte Carlo parametric bootstrap:
    %   add per-dataset Gaussian noise (sigma from fit residuals) to model
    %   predictions, refit N times, take 2.5/97.5 percentiles as 95% CI.
    fprintf('Running Monte Carlo CI estimation (%d iterations)...\n', pars.mc_n_iter);

    p_fit_tmp   = make_param_struct(x_fit, param_names);
    sigma_ds    = zeros(numel(dataset_list), 1);
    model_preds = cell(numel(dataset_list), 1);
    for d = 1:numel(dataset_list)
        ds_tmp = dataset_list{d};
        mp = ds_tmp.model_fun(p_fit_tmp, ds_tmp.indep_vars{:}, lsfcn);
        model_preds{d} = mp(:);
        r_raw = mp(:) - ds_tmp.data(:);
        sigma_ds(d) = std(r_raw);
%         sigma_ds(d) = rms(r_raw);
%         sigma_ds(d) = 1.4826 * median(abs(r_raw - median(r_raw)));
        if sigma_ds(d) == 0
            sigma_ds(d) = eps;
        end
    end

    opts_mc = optimoptions('lsqnonlin', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 16000, ...
        'MaxIterations', 1600, ...
        'UseParallel', false, ...
        'FunctionTolerance', 1e-16, ...
        'StepTolerance', 1e-16, ...
        'OptimalityTolerance', 1e-16);

    n_iter    = pars.mc_n_iter;
    n_params  = numel(param_names);
    mc_params = zeros(n_iter, n_params);

    if pars.mc_use_parallel
        parfor ii = 1:n_iter
            noisy_ds = dataset_list;
            for d = 1:numel(dataset_list)
                noisy_ds{d}.data = model_preds{d} + sigma_ds(d) * randn(size(model_preds{d}));
            end
            objfun_mc = @(x) total_residual(x, param_names, noisy_ds, lsfcn);
            mc_params(ii,:) = lsqnonlin(objfun_mc, x_fit, lb, ub, opts_mc);
        end
    else
        for ii = 1:n_iter
            noisy_ds = dataset_list;
            for d = 1:numel(dataset_list)
                noisy_ds{d}.data = model_preds{d} + sigma_ds(d) * randn(size(model_preds{d}));
            end
            objfun_mc = @(x) total_residual(x, param_names, noisy_ds, lsfcn);
            mc_params(ii,:) = lsqnonlin(objfun_mc, x_fit, lb, ub, opts_mc);
        end
    end

    ci = prctile(mc_params, [2.5 97.5], 1)';   % [n_params x 2]

else
    % Linearized asymptotic CIs via nlparci
    ci = nlparci(x_fit, residual, 'jacobian', jacobian);
end

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

    % Monte Carlo CIs for derived biexponential parameters
    ci_biexp = [];
    if use_mc
        mc_biexp = zeros(n_iter, 5);
        for ii = 1:n_iter
            ra_i  = mc_params(ii, strcmp(param_names,'ra'));
            rb_i  = mc_params(ii, strcmp(param_names,'rb'));
            mb0_i = mc_params(ii, strcmp(param_names,'mb0'));
            r_i   = mc_params(ii, strcmp(param_names,'r'));
            t2b_i = mc_params(ii, strcmp(param_names,'t2b')) * 1e-6;
            ma0_i = mc_params(ii, strcmp(param_names,'ma0'));

            R1p_i = 0.5*(ra_i + rb_i + r_i*mb0_i + r_i + sqrt((ra_i - rb_i + r_i*mb0_i - r_i)^2 ...
                + 4*r_i*mb0_i*r_i));
            R1n_i = 0.5*(ra_i + rb_i + r_i*mb0_i + r_i - sqrt((ra_i - rb_i + r_i*mb0_i - r_i)^2 ...
                + 4*r_i*mb0_i*r_i));

            mbstart_i = exp(-w1p^2 * Tp * t2b_i * adjExpFac);
            Bp_i = ((ma0_i - 1)*(ra_i - R1n_i) + (ma0_i - mbstart_i)*r_i*mb0_i)/(R1p_i - R1n_i);
            Bn_i = -((ma0_i - 1)*(ra_i - R1p_i) + (ma0_i - mbstart_i)*r_i*mb0_i)/(R1p_i - R1n_i);

            mc_biexp(ii,:) = [R1p_i, R1n_i, Bp_i, Bn_i, mbstart_i];
        end
        ci_biexp = prctile(mc_biexp, [2.5 97.5], 1)';  % [5 x 2]
    end
else
    selIRflg=false;
    ci_biexp = [];
end

%% --- Display results ---
fprintf('\n=== Fitted Parameters (95%% CI via %s) ===\n', upper(pars.ci_method));
for i = 1:numel(param_names)
    fprintf('%-10s : %.5g   (95%% CI: %.5g to %.5g)\n', ...
        param_names{i}, x_fit(i), ci(i,1), ci(i,2));
end

if selIRflg
    fprintf('\n=== Fitted Biexponential Parameters ===\n');
    for i = 1:numel(param_names_biexp)
        if ~isempty(ci_biexp)
            fprintf('%-10s : %.5g   (95%% CI: %.5g to %.5g)\n', ...
                param_names_biexp{i}, bxfit(i), ci_biexp(i,1), ci_biexp(i,2));
        else
            fprintf('%-10s : %.5g\n', ...
                param_names_biexp{i}, bxfit(i));
        end
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
        ds.plot_fun(ds.indep_vars, data, model_vals, ds.xaxislogscale);
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
        if ~isempty(ci_biexp)
            biexp_ci_lo = ci_biexp(:,1);
            biexp_ci_hi = ci_biexp(:,2);
        else
            biexp_ci_lo = NaN(5,1);
            biexp_ci_hi = NaN(5,1);
        end
        biexp_table = table(param_names_biexp', bxfit', biexp_ci_lo, biexp_ci_hi, ...
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