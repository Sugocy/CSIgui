function Tsom = getTsom(ori)
% Returns slice orientation matrix Tsom for sagital, coronal or transverse
% orientations.

switch lower(ori)
    case 'sag'
        Tsom = [0  0 -1 0; ...
                0 -1  0 0; ...
                1  0  0 1];
    case 'cor'
        Tsom = [0 -1 0 0; ...
                0  0 1 0; ...
                1  0 0 1];
    case 'tra'
        Tsom = [ 0 -1 0 0;
                -1  0 0 0;...
                 0  0 1 1];
end