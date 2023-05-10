function imgData = screenshot(area)
% Take a quick screenshot of a part of the screen using java-robot.
%
% Input: area = [left top width height]
% If no input is given, the full screen is snapshotted.

if nargin < 1, area = get(0,'screensize'); end

robot = java.awt.Robot();
pos = area; % area describes [left top width height]
rect = java.awt.Rectangle(pos(1),pos(2),pos(3),pos(4));
cap = robot.createScreenCapture(rect);

% Convert to an RGB image
rgb = typecast(cap.getRGB(0,0,cap.getWidth,cap.getHeight,[],0,cap.getWidth),'uint8');
imgData = zeros(cap.getHeight,cap.getWidth,3,'uint8');
imgData(:,:,1) = reshape(rgb(3:4:end),cap.getWidth,[])';
imgData(:,:,2) = reshape(rgb(2:4:end),cap.getWidth,[])';
imgData(:,:,3) = reshape(rgb(1:4:end),cap.getWidth,[])';