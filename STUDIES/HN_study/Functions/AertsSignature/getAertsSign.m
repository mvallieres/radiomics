function [signature] = getAertsSign(vol,mask,spatialRef,scanType)

% PARAMETERS
if strcmp(scanType,'CT')
    res = [1,1,3]; % In (Aerts,2014), all features were calculated from 1 X 1 X 3 mm^3 raw images
    levels = 0:25:3000; % Division of quantized bins in (Aerts,2014) for texture analysis
elseif strcmp(scanType,'PET') % Novelty tested in this work
    levels = 0:1:100; % Division of quantized bins for texture analysis, similarly to (Aerts,2014)
    res = [4,4,4]; % Resolution is here set up to the approximate resolution of PET scans
end

% INITIALIZATION
signature = zeros(1,4);
nLevels = numel(levels) - 1;


% STEP 0: RESAMPLING OF VOLUME AND MASK
pixelW = spatialRef.PixelExtentInWorldX;
sliceS = spatialRef.PixelExtentInWorldZ;
a = pixelW/res(1); b = pixelW/res(2); c = sliceS/res(3);
maskBox = imresize3D(mask,[],[round(double(size(mask,1))*a),round(double(size(mask,2))*b),round(double(size(mask,3))*c)],'nearest','fill');
ROIbox = imresize3D(vol,[],[round(double(size(vol,1))*a),round(double(size(vol,2))*b),round(double(size(vol,3))*c)],'cubic','fill');
ROIonly = ROIbox; ROIonly(~maskBox) = NaN;


% STEP 1: CALCULATE FEATURE 1
energy = ROIonly .^2; energy = sum(energy(~isnan(energy)));
signature(1) = energy;


% STEP 2: CALCULATE FEATURE 2
maskArea = imfill(maskBox,'holes'); maskArea = padarray(maskArea,[1,1,1],0); % Not using inner surfaces
volume = sum(maskBox(:)) * prod(res); % In mm^3
[p,q,r] = meshgrid(res(1).*(1:1:size(maskArea,2)),res(2).*(1:1:size(maskArea,1)),res(3).*(1:1:size(maskArea,3)));
[faces,vertices] = isosurface(p,q,r,maskArea,.5);
a = vertices(faces(:,2),:) - vertices(faces(:,1),:);
b = vertices(faces(:,3),:) - vertices(faces(:,1),:);
c = cross(a,b,2);
surface = 1/2 * sum(sqrt(sum(c.^2,2)));
compactness = volume/(sqrt(pi)*surface^(2/3));
signature(2) = compactness;


% STEP 3: CALCULATE FEATURE 3
ROIquant = ROIonly;
ROIquant(ROIonly < levels(1)) = 1; 
for n = 1:nLevels
    ROIquant(ROIonly>=levels(n) & ROIonly<levels(n+1)) = n;
end
ROIquant(ROIonly >= levels(end)) = nLevels;
newLevels = 1:max(ROIquant(~isnan(ROIquant)));
[GLN] = getGLN_Aerts(ROIquant,newLevels);
signature(3) = GLN;


% STEP 4: CALCULATE FEATURE 4
% - Similarly to (Aerts & Velazquez et al., 2014), MATLAB conventions are used
% - HLH directions: H1)  MATLAB = top-bottom --> axial image = anterior-posterior; L): MATLAB = left-right --> axial image = lef-right; H2): MATLAB = 3rd dimension --> axial image = inferior-superior
% - For the swt2 function, we have in DICOM XY directions: [LL,LH,HL,HH] = swt2(image,N,'wavName')
% - For the swt2 function, we have in MATLAB IJ (row,column) directions: [LL,HL,LH,HH] = swt2(image,N,'wavName')
sizeV = size(ROIbox);
if mod(sizeV(1),2)
    sizeV(1) = sizeV(1) + 1;
end
if mod(sizeV(2),2)
    sizeV(2) = sizeV(2) + 1;
end
if mod(sizeV(3),2)
    sizeV(3) = sizeV(3) + 1;
end
HLH = zeros(sizeV);
maskBox = imresize3D(maskBox,[],sizeV,'nearest','fill'); % Making sure the dimensions are even for decomposition at level 1
ROIbox = imresize3D(ROIbox,[],sizeV,'cubic','fill');
for k = 1:sizeV(3) % In this image plane (axial), we need the HL sub-band
    image = ROIbox(:,:,k);
    [~,HL,~,~] = swt2(image,1,'coif1');
    HLH(:,:,k) = HLH(:,:,k) + HL;
end
for i = 1:sizeV(1) % In this image plane (coronal), we need the HL sub-band
    image = squeeze(ROIbox(i,:,:))';
    [~,HL,~,~] = swt2(image,1,'coif1');
    for k = 1:sizeV(3)
        HLH(i,:,k) = HLH(i,:,k) + HL(k,:);
    end
end
for j = 1:sizeV(2) % In this image plane (sagittal), we need the HH sub-band
    image = squeeze(ROIbox(:,j,:));
    [~,~,~,HH] = swt2(image,1,'coif1');
    for k = 1:sizeV(3)
        HLH(:,j,k) = HLH(:,j,k) + HH(:,k);
    end
end
HLH = HLH/3;

ROIquant = HLH; ROIquant(~maskBox) = NaN; 
Ng = numel(newLevels); [ROIquant,levels] = uniformQuantization(ROIquant,Ng);
[GLN] = getGLN_Aerts(ROIquant,levels);
signature(4) = GLN;

end