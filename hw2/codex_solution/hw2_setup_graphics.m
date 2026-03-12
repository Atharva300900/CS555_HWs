function hw2_setup_graphics()
% Keep plotting non-interactive in Octave/MATLAB.
setenv('GNUTERM', 'pngcairo');
set(0, 'defaultfigurevisible', 'off');
end
