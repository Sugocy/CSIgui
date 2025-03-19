function moveMouse(pos)

scrsz = get(0, 'screensize'); w = scrsz(3); h = scrsz(4);

import java.awt.Robot;
import java.awt.event.*;

% There is a problem with mouseMove wrt MATLAB/JAVA/Win and high-dpi
% monitors.Matlab allows mouse-movent via the set-Pointerlocation.
% 
% 1. Try Java-way, check location.
mouse = Robot; mouse.mouseMove(pos(1), h-pos(2));
% 2. Check location
newPos = get(0,'PointerLocation');
if newPos(2) ~= (pos(2) - mod(pos(2),1)    )
    % 3. Not correct location - set via Matlab.
    set(0,'PointerLocation',[pos(1) pos(2)]);
end

% Wait a ms.
pause(10e-3);

% Click
mouse.mousePress(InputEvent.BUTTON1_MASK);
mouse.mouseRelease(InputEvent.BUTTON1_MASK);