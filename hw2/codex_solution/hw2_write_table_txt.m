function hw2_write_table_txt(filename, headers, rows)
% Write a simple tab-separated text table.
fid = fopen(filename, 'w');
if fid < 0
    error('Failed to open %s for writing.', filename);
end
cleanup = onCleanup(@() fclose(fid));

for j = 1:numel(headers)
    fprintf(fid, '%s', headers{j});
    if j < numel(headers)
        fprintf(fid, '\t');
    else
        fprintf(fid, '\n');
    end
end

for i = 1:size(rows, 1)
    for j = 1:size(rows, 2)
        fprintf(fid, '%s', rows{i, j});
        if j < size(rows, 2)
            fprintf(fid, '\t');
        else
            fprintf(fid, '\n');
        end
    end
end
end
