function about_text = loadAbout()
% Opens the about text from file for CSIgui.

% Get file ID
fid = fopen('about.txt','r');

n = 1; % Iterate
% Read each line.
about_text = [];
while ~feof(fid)
    line = fgetl(fid);              % Get line
    about_text = [about_text line]; % Add line to text
    n = n + 1;                      % Iteration up
end
fclose(fid); % Close file.
