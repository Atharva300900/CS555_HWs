function paths = hw2_paths()
% Output locations for figures, tables, and cached data.
root = fileparts(mfilename('fullpath'));
paths.root = root;
paths.output = fullfile(root, 'output');
paths.figures = fullfile(paths.output, 'figures');
paths.tables = fullfile(paths.output, 'tables');
paths.data = fullfile(paths.output, 'data');

ensure_dir(paths.output);
ensure_dir(paths.figures);
ensure_dir(paths.tables);
ensure_dir(paths.data);
end

function ensure_dir(path_name)
if exist(path_name, 'dir') ~= 7
    mkdir(path_name);
end
end
