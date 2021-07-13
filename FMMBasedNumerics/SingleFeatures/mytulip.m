xvals=0.25:.01:3;
yvals=(0.25:.01:3).^2+(0.25:.01:3).^4;
xvals=[xvals 1.5 0 -1.5];
yvals=[yvals 70 3.^2+3.^4 70];
xvals=[xvals -3:.01:-0.25];
yvals=[yvals (-3:.01:-0.25).^2+(-3:.01:-0.25).^4];
xvals=[xvals -0.25 0.25 0.25];
yvals=[yvals -(3.^2+3.^4) -(3.^2+3.^4) 0];
xvals=xvals/max(xvals(:))/2;
yvals=yvals/max(yvals(:));
%close all; plot(xvals,yvals); axis equal