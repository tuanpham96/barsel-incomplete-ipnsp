clc; clear; %close all; 
v = 0:0.01:1;
[x,y] = meshgrid(v);
dw = y.*x - 1/5*(x.^2 + y.^2); 

cmap = return_colorbrewer('RdBu_r', 100);
figure; hold on; colormap(cmap); 
image(v, v, dw', 'AlphaData', ~isnan(dw'));
daspect([1,1,1]); 
caxis([-1,1]*0.2);xlim([0,1]);ylim([0,1]);
xlabel('x'); ylabel('y');
colorbar;
despline;
