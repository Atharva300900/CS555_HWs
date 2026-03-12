function all_results = run_all()
% Run HW2 drivers a-e.
all_results.q1a = solve_q1a();
all_results.q1b = solve_q1b();
all_results.q1c = solve_q1c();
all_results.q1d = solve_q1d();
all_results.q1e = solve_q1e();
save(fullfile(hw2_paths().data, 'run_all_results.mat'), 'all_results');
end
