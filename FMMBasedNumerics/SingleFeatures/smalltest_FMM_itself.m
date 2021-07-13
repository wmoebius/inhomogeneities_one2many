clear all;  %#ok<CLALL>
close all;
clc;

% clear folder
delete('TEMP/*')

% set up features and speedmatrix
%[X,Y]=meshgrid(-1:0.01:1,-1:0.01:1);
%speedmatrix=1.0+X.^2+Y.^2;

myellipse2
[X,Y]=meshgrid(-1.5:0.01:3,-2.5:0.02:2.5);
myinpol=inpolygon(X(:),Y(:),yvals,xvals);
speedmatrix=1.0+0.5*reshape(myinpol,size(X,1),size(X,2));

% FMM
dlmwrite('TEMP/myX.txt',X,' ');
dlmwrite('TEMP/myY.txt',Y,' ');
dlmwrite('TEMP/myspeedmatrix.txt',speedmatrix,' ');
system('source ~/.bash_profile; source activate skfmm; python singlefeaturefmm.py');

subplot(2,2,1)
speedmatrix=dlmread('TEMP/speedmatrix.dat');
imagesc(1./speedmatrix)
axis equal
colorbar
title('input: inverse speed matrix')

subplot(2,2,2)
tmatrix=dlmread('TEMP/nanfilledtm.dat');
imagesc(tmatrix)
axis equal
colorbar
title('ouput skfmm: arrival time')

subplot(2,2,3)
[fx,fy]=gradient(tmatrix);
imagesc(sqrt(fx.^2+fy.^2));
axis equal
colorbar
title('gradient of arrival time (ouput skfmm)')

subplot(2,2,4)
[fx,fy]=gradient(tmatrix);
imagesc(sqrt(fx.^2+fy.^2).*speedmatrix*100);
axis equal
colorbar
title('gradient of arrival time (ouput skfmm) * speedmatrix * 100')

figure

subplot(1,3,1)
tmatrix=dlmread('TEMP/nanfilledtm.dat');
imagesc(tmatrix)
axis equal
colorbar
title('ouput skfmm: arrival time')

subplot(1,3,2)
[~,matlabfmmoutput]=imsegfmm(speedmatrix,X==min(X(:)),.1);
imagesc(matlabfmmoutput)
axis equal
colorbar
title('ouput Matlab: arrival time')

subplot(1,3,3)
imagesc(tmatrix/max(tmatrix(:))-matlabfmmoutput)
axis equal
colorbar
title('differece skfmm Matlab')