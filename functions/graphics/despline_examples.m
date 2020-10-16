%% Example 1: with both log scale 
figure; hold on; 
plot((10:100).^2,10:100); 
xlim([10, 0.5e5]); ylim([1, 200]); 
set(gca, 'yscale', 'log'); set(gca, 'xscale', 'log')
despline(0.5);

%% Example 2: only 1 log 
figure; hold on; 
plot((10:100).^2,10:100); 
ylim([1, 200]); 
set(gca, 'yscale', 'log'); 
despline(1);

%% Example 3: all axes of a figure 
figure; 
for i = 1:10; subplot(3,4,i); end
despline(gcf,2)

%% Example 4: selective 
figure; 
a = subplot(312); 
b = subplot(311); c = subplot(313); 
despline([a,b],[1,2])