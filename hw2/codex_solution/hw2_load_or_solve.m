function results = hw2_load_or_solve(N, grid_type, form, target_dt, t_save)
% Cache solves so the drivers can reuse common cases.
if nargin < 5
    t_save = [];
end
cfg = hw2_config();
paths = hw2_paths();
[solver_cfl, target_dt] = hw2_effective_cfl(N, grid_type, target_dt);

if isempty(t_save)
    save_tag = 'nosave';
else
    save_tag = sprintf('save%d', numel(t_save));
end

cache_version = 'v2';
dt_str = sprintf('%.10e', target_dt);
dt_str = strrep(dt_str, '.', 'p');
dt_str = strrep(dt_str, '+', '');
dt_str = strrep(dt_str, '-', 'm');
cache_name = sprintf('%s_N%d_%s_%s_dt%s_%s.mat', cache_version, N, lower(grid_type), lower(form), dt_str, save_tag);
cache_path = fullfile(paths.data, cache_name);

if exist(cache_path, 'file') == 2
    tmp = load(cache_path);
    results = tmp.results;
    return;
end

results = burgers_solve(N, cfg.nu, cfg.T, grid_type, form, solver_cfl, t_save);
results.grid_type = grid_type;
results.form = form;
results.target_dt = target_dt;
results.solver_cfl = solver_cfl;
save(cache_path, 'results');
end
