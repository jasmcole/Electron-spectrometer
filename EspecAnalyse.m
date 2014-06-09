function output = EspecAnalyse(month, day, run, shot)

clf

I1 = GetFilename(month, day, run, shot, 'Espec1.tif');
I1 = double(imread(I1));
I2 = GetFilename(month, day, run, shot, 'Espec2.tif');
I2 = double(imread(I2));

if(month == xx && day == xx)
    q1 = load('xx_Espec1_warping.mat');
    q1 = q1.q;
    q2 = load('xx_Espec2_warping.mat');
    q2 = q2.q;
    voffset = 46;
    I1(1:end-voffset,:) = I1(voffset+1:end,:);
end

subplot(2,3,1)
imagesc(q1.u, q1.v, I1); axis image xy; caxis([0 5000])
hold on
plot([q1.screenu' q1.screenu'], [q1.screenv' q1.screenv'], 'w')
hold off
title('Screen 1 - Original Image', 'FontSize', 14)
set(gca, 'XTick', []);
set(gca, 'YTick', []);

subplot(2,3,2)
imagesc(q2.u, q2.v, I2); axis image xy; caxis([0 5000])
hold on
plot([q2.screenu' q2.screenu'], [q2.screenv' q2.screenv'], 'w')
hold off
title('Screen 2 - Original Image', 'FontSize', 14)
set(gca, 'XTick', []);
set(gca, 'YTick', []);
drawnow

%Transform image, preserving counts overall sum(sum(I)) is conserved
I1warp = ProjectiveTransformImage(I1,q1);
I2warp = ProjectiveTransformImage(I2,q2);

%Find middle 80% of lanex axis
screenheight1 = max(q1.screeny);
screenheight2 = max(q2.screeny);
[~, topsection1] = min(abs(q1.y - (9/10)*screenheight1));
[~, botsection1] = min(abs(q1.y - (1/10)*screenheight1));
[~, topsection2] = min(abs(q2.y - (9/10)*screenheight2));
[~, botsection2] = min(abs(q2.y - (1/10)*screenheight2));

[~, screentop1] = min(abs(q1.y - screenheight1));
[~, screentop2] = min(abs(q2.y - screenheight2));

%Subtract background
[~, zero1]  = min(abs(q1.y));
[~, zero2]  = min(abs(q2.y));

bg1 = 0.5*(mean(I1warp(zero1:botsection1,:)) + mean(I1warp(topsection1:screentop1,:)));
bg2 = 0.5*(mean(I2warp(zero2:botsection2,:)) + mean(I2warp(topsection2:screentop2,:)));
bgspline1 = csaps(1:length(bg1), bg1, 0.1, 1:length(bg1));
bgspline2 = csaps(1:length(bg2), bg2, 0.1, 1:length(bg2));
bg1 = ones(length(I1warp(:,1)),1)*bgspline1;
bg2 = ones(length(I2warp(:,1)),1)*bgspline2;
I1warp = I1warp - bg1;
I2warp = I2warp - bg2;

%Plot initial and warped images
subplot(2,3,4)
imagesc(q1.x/q1.pixpermm, q1.y/q1.pixpermm, I1warp); axis image xy; caxis([0 2000])
hold on
plot([q1.screenx/q1.pixpermm q1.screenx/q1.pixpermm], [q1.screeny/q1.pixpermm q1.screeny/q1.pixpermm], 'w')
line([0 max(q1.screenx/q1.pixpermm)], [q1.y(botsection1)/q1.pixpermm q1.y(botsection1)/q1.pixpermm], 'color', 'white', 'linestyle', '--')
line([0 max(q1.screenx/q1.pixpermm)], [q1.y(topsection1)/q1.pixpermm q1.y(topsection1)/q1.pixpermm], 'color', 'white', 'linestyle', '--')
hold off
title('Screen 1 - Warped Image', 'FontSize', 14)
xlabel('x /mm', 'FontSize', 14)
ylabel('y /mm', 'FontSize', 14)

subplot(2,3,5)
imagesc(q2.x/q2.pixpermm, q2.y/q2.pixpermm, I2warp); axis image xy; caxis([0 2000])
hold on
plot([q2.screenx/q2.pixpermm q2.screenx/q2.pixpermm], [q2.screeny/q2.pixpermm q2.screeny/q2.pixpermm], 'w')
line([0 max(q2.screenx/q2.pixpermm)], [q2.y(botsection2)/q2.pixpermm q2.y(botsection2)/q2.pixpermm], 'color', 'white', 'linestyle', '--')
line([0 max(q2.screenx/q2.pixpermm)], [q2.y(topsection2)/q2.pixpermm q2.y(topsection2)/q2.pixpermm], 'color', 'white', 'linestyle', '--')
hold off
title('Screen 2 - Warped Image', 'FontSize', 14)
xlabel('x /mm', 'FontSize', 14)
ylabel('y /mm', 'FontSize', 14)




P1 = 1:length(I1warp(1,:));
P2 = 1:length(I2warp(1,:));

%d is distance from high energy edge
d1 = double(q1.x)/q1.pixpermm;
d2 = double(q2.x)/q2.pixpermm;
[~,screen1right] = min(abs(q1.x - max(q1.screenx)));
[~,screen2right] = min(abs(q2.x - max(q2.screenx)));
d1 = q1.x(screen1right)/q1.pixpermm - d1;
d2 = q2.x(screen2right)/q2.pixpermm - d2;

% Magnet tracking might need to be re-done with better map
load('ScreenCalibration.mat')
E_MeV1 = csaps(daxis1, Eaxis1, 0.9, d1);
Ewidth1 = gradient(E_MeV1);

E_MeV2 = csaps(daxis2, Eaxis2, 0.9, d2);
Ewidth2 = gradient(E_MeV2);

counts1 = sum(I1warp(botsection1:topsection1,:));
counts2 = sum(I2warp(botsection2:topsection2,:));
totalCounts1 = sum(counts1);
totalCounts2 = sum(counts2);
totalCharge1 = totalCounts1 * 2.89e-5;
totalCharge2 = totalCounts2 * 2.89e-5;
disp(['Total Charge on screen 1 ' num2str(totalCharge1) 'pC'])
disp(['Total Charge on screen 2 ' num2str(totalCharge2) 'pC'])
countsperMeV1 = counts1./Ewidth1;
countsperMeV2 = counts2./Ewidth2;

%Find edges of lanex screen
[~, high_E_edge_index1] = min(abs(q1.x - max(q1.screenx)));
[~, low_E_edge_index1] =  min(abs(q1.x));
[~, high_E_edge_index2] = min(abs(q2.x - max(q2.screenx)));
[~, low_E_edge_index2] =  min(abs(q2.x));

countsperMeV1 = countsperMeV1(low_E_edge_index1:high_E_edge_index1);
E_MeV1 = E_MeV1(low_E_edge_index1:high_E_edge_index1);
countsperMeV2 = countsperMeV2(low_E_edge_index2:high_E_edge_index2);
E_MeV2 = E_MeV2(low_E_edge_index2:high_E_edge_index2);

countsperMeV1 = csaps(E_MeV1, countsperMeV1, 0.6, E_MeV1);
countsperMeV2 = csaps(E_MeV2, countsperMeV2, 0.6, E_MeV2);

%pC and MeV
Qaxis1 = countsperMeV1 * totalCharge1 / trapz(E_MeV1, countsperMeV1);
Eaxis1 = E_MeV1;
Qaxis2 = countsperMeV2 * totalCharge2 / trapz(E_MeV2, countsperMeV2);
Eaxis2 = E_MeV2;

subplot(2,3,[3 6])
plot(Eaxis1, Qaxis1, Eaxis2, Qaxis2)
xlim([100 1100])
ylabel('Charge / pCMeV^{-1}', 'FontSize', 14)
xlabel('Energy /MeV', 'FontSize', 14)
title('Spectrum', 'FontSize', 14)
legend('Screen 1', 'Screen 2')
set(gcf, 'Color', 'white')

mTextBox = uicontrol('style','text');
set(mTextBox,'String',[num2str(day) '/' num2str(month) '/2012 Run ' num2str(run) ', shot ' num2str(shot)])
set(mTextBox, 'Position', [100 0 300 30],'BackgroundColor', 'white','HorizontalAlignment', 'left', 'FontSize', 14)

output.Eaxis1 = Eaxis1;
output.Eaxis2 = Eaxis2;
output.Qaxis1 = Qaxis1;
output.Qaxis2 = Qaxis2;
output.screen1 = I1warp;
output.screen2 = I2warp;

end