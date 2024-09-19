function moveMouse(pos)

scrsz = get(0, 'screensize'); w = scrsz(3); h = scrsz(4);

import java.awt.Robot;
import java.awt.event.*;

mouse = Robot;
mouse.mouseMove(pos(1), h-pos(2));

pause(0.0001);

mouse.mousePress(InputEvent.BUTTON1_MASK);
mouse.mouseRelease(InputEvent.BUTTON1_MASK);