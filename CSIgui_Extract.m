% Read file
txt = fileread('CSIgui.m');
% Split lines
newStr = splitlines(txt);
% Find functions
a = contains(newStr,'function ');
func_total = sum(a==1); func_lines = newStr(a==1);
% Remove possible comments with "function" 
b = contains(func_lines, '%'); func_lines = func_lines(b == 0);


fid = fopen('CSIgui_FunctionsList.txt','wt+');
for kk = 1:size(func_lines)
   fprintf(fid, [func_lines{kk} '\n']);
    
end
fclose(fid);