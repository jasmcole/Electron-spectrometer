function q = GetProjectiveTransformImage

clf
[FileName,PathName,FilterIndex] = uigetfile('*.*');
inputImage = [PathName FileName];

imageSize = size(inputImage);

if(imageSize(1) == 1)
    if(strcmp(inputImage(end-2:end), 'raw'))
        inputImage = ReadRAW16bit(inputImage, 1280, 960);
    else
        inputImage = imread(inputImage);
    end
end

inputImage = double(inputImage);
figure(1)
imagesc(inputImage)
axis image xy
colormap(gray)
xlabel('u')
ylabel('v')
title('Input Image')
imageSize = size(inputImage);
udata = [1 imageSize(2)];
vdata = [1 imageSize(1)];
disp('Start at the bottom left and go clockwise')
disp('Double-click inside when done')
title('Start at the bottom left and go clockwise. Double-click inside when done')
imcontrast
figure(1)
hpoly = impoly(gca);
pos = wait(hpoly);
x = pos(:,1);
y = pos(:,2);

lwidth  = inputdlg('Lanex width (mm)');
lheight = inputdlg('Lanex height (mm)');
pixpermm = inputdlg('Pixels per mm in output image');

lwidth = str2num(pixpermm{1,1})*str2num(lwidth{1,1});
lheight = str2num(pixpermm{1,1})*str2num(lheight{1,1});
pixpermm = str2num(pixpermm{1,1});

xp = [0 0 lwidth lwidth];
yp = [0 lheight lheight 0];

subplot(2,2,1)
imagesc(inputImage)

tform = maketform('projective',[ x(1) y(1); x(2) y(2); x(3) y(3); x(4) y(4)],...
                               [ xp(1) yp(1); xp(2) yp(2); xp(3) yp(3); xp(4) yp(4)]);
                           
u = 1:1:imageSize(2);
v = 1:1:imageSize(1);
for n = 1:length(u)
    for m = 1:length(v)
        J(m,n) = ProjectiveTransformJacobian(u(n), v(m), tform);
    end
end                           
                           
[outputImage,xdata,ydata] = imtransform(inputImage, tform, 'bicubic', ...
                                        'udata', udata,...
                                        'vdata', vdata,...
                                        'XYScale', 1);
[Jnew,xdata,ydata]        = imtransform(J, tform, 'bicubic', ...
                                        'udata', udata,...
                                        'vdata', vdata,...
                                        'XYScale', 1);
                                    
xaxis = round(xdata(1)):1:round(xdata(1))+length(outputImage(1,:))-1;
yaxis = round(ydata(1)):1:round(ydata(1))+length(outputImage(:,1))-1;

outputImage = outputImage./Jnew;
outputImage(isnan(outputImage)) = 0;
                                    
subplot(2,2,2)
imagesc(xaxis,yaxis,outputImage);
hold on
plot([xp xp], [yp yp], 'w')
hold off
axis image xy
xlabel('x')
ylabel('y')
title('Output Image')

subplot(2,2,3)
imagesc(J); axis image xy
xlabel('u')
ylabel('v')
title('Jacobian (u,v)')
subplot(2,2,4)
imagesc(xaxis, yaxis, Jnew); axis image xy
xlabel('u')
ylabel('v')
title('Jacobian (x,y)')

%q.I = inputImage;
%q.I2 = outputImage;
%q.J = J;
q.J2 = Jnew;
q.tform = tform;
q.x = xaxis;
q.y = yaxis;
q.u = u;
q.v = v;
q.screenu = x;
q.screenv = y;
q.screenx = xp;
q.screeny = yp;
q.pixpermm = pixpermm;

[FileName,PathName,FilterIndex] = uiputfile('*.mat', 'Save calibration as ');
outputFile = [PathName FileName];
save(outputFile, 'q')

end