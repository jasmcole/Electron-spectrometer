function outputImage = ProjectiveTransformImage(inputImage, q)

imageSize = size(inputImage);

if(imageSize(1) == 1)
    inputImage = imread(inputImage);
end

inputImage = double(inputImage);

imageSize = size(inputImage);
udata = [1 imageSize(2)];
vdata = [1 imageSize(1)];

[outputImage,xdata,ydata] = imtransform(inputImage, q.tform, 'bicubic', ...
                                        'udata', udata,...
                                        'vdata', vdata,...
                                        'XYScale', 1);
outputImage = outputImage./q.J2;
outputImage(isnan(outputImage)) = 0;

% subplot(1,2,1)
% set(gca, 'FontSize', 14)
% imagesc(q.u,q.v,inputImage);axis image xy
% hold on
% plot([q.screenu' q.screenu'], [q.screenv' q.screenv'], 'w')
% hold off
% subplot(1,2,2)
% set(gca, 'FontSize', 14)
% imagesc(q.x, q.y, outputImage); axis image xy
% hold on
% plot([q.screenx q.screenx], [q.screeny q.screeny], 'w')
% hold off
% axis image xy

end