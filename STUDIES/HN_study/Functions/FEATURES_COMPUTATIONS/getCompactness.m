function [compactness] = getCompactness(ROIonly,pixelW,sliceS)

% CITE Aerts AS CO-AUTHOR HERE.


mask = ~isnan(ROIonly); % Find mask covering the ROI

% ISOTROPIC RESAMPLING
sFactor = sliceS/pixelW; % scaling factor
mask = imresize3D(mask,[],[round(double(size(mask,1))),round(double(size(mask,2))),round(double(size(mask,3))*sFactor)],'nearest','fill'); % Now isotropically resampled to pixelW

% COMPACTNESS COMPUTATION
maskArea = imfill(mask,'holes'); maskArea = padarray(maskArea,[1,1,1],0); % Not using inner surfaces
volume = sum(mask(:)) * pixelW * pixelW * pixelW; % In mm^3
[p,q,r] = meshgrid(pixelW.*(1:1:size(maskArea,2)),pixelW.*(1:1:size(maskArea,1)),pixelW.*(1:1:size(maskArea,3)));
[faces,vertices] = isosurface(p,q,r,maskArea,.5);
a = vertices(faces(:,2),:) - vertices(faces(:,1),:);
b = vertices(faces(:,3),:) - vertices(faces(:,1),:);
c = cross(a,b,2);
surface = 1/2 * sum(sqrt(sum(c.^2,2)));
compactness = volume/(sqrt(pi)*surface^(3/2)); % Dimensionless version, different from (Aerts et al., 2014)

end