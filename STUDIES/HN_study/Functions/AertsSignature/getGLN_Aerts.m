function [GLN] = getGLN_Aerts(ROIonly,levels)
% Same method as Aerts, 2014


% PRELIMINARY
nLevel = length(levels);
levelTemp = max(levels)+1;
ROIonly(isnan(ROIonly)) = levelTemp; % Last row needs to be taken out of the GLRLM
levels = [levels,levelTemp]; NL = length(levels) - 1;
sizeV = size(ROIonly);
numInit = ceil(max(sizeV)*sqrt(3)); % Max run length
GLN = 0;


%%%%%%%% START COMPUTATION %%%%%%%%

% Directions [1,0,0], [0 1 0], [1 1 0] and [-1 1 0] : 2D directions
% (x:right-left, y:top-bottom, z:3rd dimension)  
nComp = sizeV(3); % We can add-up the GLRLMs taken separately in every image in the x-y plane

% [1,0,0]
GLRLM = zeros(NL+1,numInit);
for i = 1:nComp
    image = ROIonly(:,:,i);
    uniqueIm = unique(image);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image;
    for j = 1:NLtemp
        indexRow(j) = find(uniqueIm(j)==levels);
        image(temp==uniqueIm(j)) = j;
    end
    GLRLMtemp = rle_0(image,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN

% [0 1 0]
GLRLM = zeros(NL+1,numInit);
for i = 1:nComp
    image = ROIonly(:,:,i);
    uniqueIm = unique(image);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image;
    for j = 1:NLtemp
        indexRow(j) = find(uniqueIm(j)==levels);
        image(temp==uniqueIm(j)) = j;
    end
    GLRLMtemp = rle_0(image',NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN

% [1 1 0]
GLRLM = zeros(NL+1,numInit);
for i = 1:nComp
    image = ROIonly(:,:,i);
    uniqueIm = unique(image);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image;
    for j = 1:NLtemp
        indexRow(j) = find(uniqueIm(j)==levels);
        image(temp==uniqueIm(j)) = j;
    end
    seq = zigzag(image);
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN

% [-1 1 0]
GLRLM = zeros(NL+1,numInit);
for i = 1:nComp
    image = ROIonly(:,:,i);
    uniqueIm = unique(image);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image;
    for j = 1:NLtemp
        indexRow(j) = find(uniqueIm(j)==levels);
        image(temp==uniqueIm(j)) = j;
    end
    seq = zigzag(fliplr(image));
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN


% Directions [0,0,1], [1 0 1] and [-1 0 1]
% (x:right-left, y:top-bottom, z:3rd dimension)
nComp = sizeV(1); % We can add-up the GLRLMs taken separately in every image in the x-z plane

% [0,0,1]
GLRLM = zeros(NL+1,numInit);
image = zeros(sizeV(3),sizeV(2));
for i = 1:nComp
    for j = 1:sizeV(3)
        image(j,1:end) = ROIonly(i,1:end,j);
    end
    uniqueIm = unique(image);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image;
    for j=1:NLtemp
        indexRow(j) = find(uniqueIm(j)==levels);
        image(temp==uniqueIm(j)) = j;
    end
    GLRLMtemp = rle_0(image',NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN

% [1 0 1]
GLRLM = zeros(NL+1,numInit);
image = zeros(sizeV(3),sizeV(2));
for i = 1:nComp
    for j = 1:sizeV(3)
        image(j,1:end) = ROIonly(i,1:end,j);
    end
    uniqueIm = unique(image);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image;
    for j=1:NLtemp
        indexRow(j) = find(uniqueIm(j)==levels);
        image(temp==uniqueIm(j)) = j;
    end
    seq = zigzag(image);
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN

% [-1 0 1]
GLRLM = zeros(NL+1,numInit);
image = zeros(sizeV(3),sizeV(2));
for i = 1:nComp
    for j = 1:sizeV(3)
        image(j,1:end) = ROIonly(i,1:end,j);
    end
    uniqueIm = unique(image);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image;
    for j=1:NLtemp
        indexRow(j) = find(uniqueIm(j)==levels);
        image(temp==uniqueIm(j)) = j;
    end
    seq = zigzag(fliplr(image));
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN


% Directions [0,1,1] and [0 -1 1]
% (x:right-left, y:top-bottom, z:3rd dimension)
nComp = sizeV(2); % We can add-up the GLRLMs taken separately in every image in the y-z plane

% [0,1,1]
GLRLM = zeros(NL+1,numInit);
image = zeros(sizeV(1),sizeV(3));
for i = 1:nComp
    for j = 1:sizeV(3)
        image(1:end,j) = ROIonly(1:end,i,j);
    end
    uniqueIm = unique(image);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image;
    for j=1:NLtemp
        indexRow(j) = find(uniqueIm(j)==levels);
        image(temp==uniqueIm(j)) = j;
    end
    seq = zigzag(image);
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN

% [0 -1 1]
GLRLM = zeros(NL+1,numInit);
image = zeros(sizeV(1),sizeV(3));
for i = 1:nComp
    for j = 1:sizeV(3)
        image(1:end,j) = ROIonly(1:end,i,j);
    end
    uniqueIm = unique(image);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image;
    for j=1:NLtemp
        indexRow(j) = find(uniqueIm(j)==levels);
        image(temp==uniqueIm(j)) = j;
    end
    seq = zigzag(fliplr(image));
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN


% Four corners: [1,1,1], [-1,1,1], [-1,1,-1], [1,1,-1]
% (x:right-left, y:top-bottom, z:3rd dimension)
GLRLM = zeros(NL+1,numInit);
image = zeros(sizeV(3),sizeV(2));
temp = rand(sizeV(3),sizeV(2));
diagTemp = spdiags(temp);
szDiag = size(diagTemp);
diagMat1 = zeros(szDiag(1),szDiag(2),sizeV(1));
diagMat2 = zeros(size(diagTemp,1),size(diagTemp,2),sizeV(1));
for i = 1:sizeV(1)
    for j = 1:sizeV(3)
        image(j,1:end) = ROIonly(i,1:end,j);
    end
    try
        diagMat1(:,:,i)=spdiags(image);
    catch
        % Add a column at the beginning to prevent errors
        temp=spdiags(image);
        numberDiff=abs(size(temp,2)-size(diagMat1,2));
        if mod(numberDiff,2) % Odd difference number
            temp=padarray(temp,[0,(numberDiff+1)/2,0],0);
            diagMat1(:,:,i)=temp(:,1:end-1);
        else
            diagMat1(:,:,i)=padarray(temp,[0,numberDiff/2,0],0);
        end
    end
    try
        diagMat2(:,:,i)=spdiags(fliplr(image));
    catch
        % Add a column at the beginning to prevent errors
        temp = spdiags(fliplr(image));
        numberDiff = abs(size(temp,2)-size(diagMat2,2));
        if mod(numberDiff,2) % Odd difference number
            temp = padarray(temp,[0,(numberDiff+1)/2,0],0);
            diagMat2(:,:,i) = temp(:,1:end-1);
        else
            diagMat2(:,:,i) = padarray(temp,[0,numberDiff/2,0],0);
        end
    end
end
for j = 1:szDiag(2)
    index = (diagMat1(:,j,1)~=0);
    nTemp = sum(index);
    image1 = zeros(sizeV(1),nTemp);
    image2 = zeros(sizeV(1),nTemp);
    for k = 1:sizeV(1)
        image1(k,1:nTemp) = diagMat1(index(1:end),j,k)';
        image2(k,1:nTemp) = diagMat2(index(1:end),j,k)';
    end
    uniqueIm = unique(image1);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image1;
    for i = 1:NLtemp
        indexRow(i) = find(uniqueIm(i)==levels);
        image1(temp==uniqueIm(i)) = i;
    end
    seq = zigzag(image1);
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN

GLRLM = zeros(NL+1,numInit);
image = zeros(sizeV(3),sizeV(2));
temp = rand(sizeV(3),sizeV(2));
diagTemp = spdiags(temp);
szDiag = size(diagTemp);
diagMat1 = zeros(szDiag(1),szDiag(2),sizeV(1));
diagMat2 = zeros(size(diagTemp,1),size(diagTemp,2),sizeV(1));
for i = 1:sizeV(1)
    for j = 1:sizeV(3)
        image(j,1:end) = ROIonly(i,1:end,j);
    end
    try
        diagMat1(:,:,i)=spdiags(image);
    catch
        % Add a column at the beginning to prevent errors
        temp=spdiags(image);
        numberDiff=abs(size(temp,2)-size(diagMat1,2));
        if mod(numberDiff,2) % Odd difference number
            temp=padarray(temp,[0,(numberDiff+1)/2,0],0);
            diagMat1(:,:,i)=temp(:,1:end-1);
        else
            diagMat1(:,:,i)=padarray(temp,[0,numberDiff/2,0],0);
        end
    end
    try
        diagMat2(:,:,i)=spdiags(fliplr(image));
    catch
        % Add a column at the beginning to prevent errors
        temp = spdiags(fliplr(image));
        numberDiff = abs(size(temp,2)-size(diagMat2,2));
        if mod(numberDiff,2) % Odd difference number
            temp = padarray(temp,[0,(numberDiff+1)/2,0],0);
            diagMat2(:,:,i) = temp(:,1:end-1);
        else
            diagMat2(:,:,i) = padarray(temp,[0,numberDiff/2,0],0);
        end
    end
end
for j = 1:szDiag(2)
    index = (diagMat1(:,j,1)~=0);
    nTemp = sum(index);
    image1 = zeros(sizeV(1),nTemp);
    image2 = zeros(sizeV(1),nTemp);
    for k = 1:sizeV(1)
        image1(k,1:nTemp) = diagMat1(index(1:end),j,k)';
        image2(k,1:nTemp) = diagMat2(index(1:end),j,k)';
    end
    uniqueIm = unique(image1);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image1;
    for i = 1:NLtemp
        indexRow(i) = find(uniqueIm(i)==levels);
        image1(temp==uniqueIm(i)) = i;
    end
    seq = zigzag(fliplr(image1));
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN

GLRLM = zeros(NL+1,numInit);
image = zeros(sizeV(3),sizeV(2));
temp = rand(sizeV(3),sizeV(2));
diagTemp = spdiags(temp);
szDiag = size(diagTemp);
diagMat1 = zeros(szDiag(1),szDiag(2),sizeV(1));
diagMat2 = zeros(size(diagTemp,1),size(diagTemp,2),sizeV(1));
for i = 1:sizeV(1)
    for j = 1:sizeV(3)
        image(j,1:end) = ROIonly(i,1:end,j);
    end
    try
        diagMat1(:,:,i)=spdiags(image);
    catch
        % Add a column at the beginning to prevent errors
        temp=spdiags(image);
        numberDiff=abs(size(temp,2)-size(diagMat1,2));
        if mod(numberDiff,2) % Odd difference number
            temp=padarray(temp,[0,(numberDiff+1)/2,0],0);
            diagMat1(:,:,i)=temp(:,1:end-1);
        else
            diagMat1(:,:,i)=padarray(temp,[0,numberDiff/2,0],0);
        end
    end
    try
        diagMat2(:,:,i)=spdiags(fliplr(image));
    catch
        % Add a column at the beginning to prevent errors
        temp = spdiags(fliplr(image));
        numberDiff = abs(size(temp,2)-size(diagMat2,2));
        if mod(numberDiff,2) % Odd difference number
            temp = padarray(temp,[0,(numberDiff+1)/2,0],0);
            diagMat2(:,:,i) = temp(:,1:end-1);
        else
            diagMat2(:,:,i) = padarray(temp,[0,numberDiff/2,0],0);
        end
    end
end
for j = 1:szDiag(2)
    index = (diagMat1(:,j,1)~=0);
    nTemp = sum(index);
    image1 = zeros(sizeV(1),nTemp);
    image2 = zeros(sizeV(1),nTemp);
    for k = 1:sizeV(1)
        image1(k,1:nTemp) = diagMat1(index(1:end),j,k)';
        image2(k,1:nTemp) = diagMat2(index(1:end),j,k)';
    end
    uniqueIm = unique(image2);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image2;
    for i = 1:NLtemp
        indexRow(i) = find(uniqueIm(i)==levels);
        image2(temp==uniqueIm(i)) = i;
    end
    seq = zigzag(image2);
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN

GLRLM = zeros(NL+1,numInit);
image = zeros(sizeV(3),sizeV(2));
temp = rand(sizeV(3),sizeV(2));
diagTemp = spdiags(temp);
szDiag = size(diagTemp);
diagMat1 = zeros(szDiag(1),szDiag(2),sizeV(1));
diagMat2 = zeros(size(diagTemp,1),size(diagTemp,2),sizeV(1));
for i = 1:sizeV(1)
    for j = 1:sizeV(3)
        image(j,1:end) = ROIonly(i,1:end,j);
    end
    try
        diagMat1(:,:,i)=spdiags(image);
    catch
        % Add a column at the beginning to prevent errors
        temp=spdiags(image);
        numberDiff=abs(size(temp,2)-size(diagMat1,2));
        if mod(numberDiff,2) % Odd difference number
            temp=padarray(temp,[0,(numberDiff+1)/2,0],0);
            diagMat1(:,:,i)=temp(:,1:end-1);
        else
            diagMat1(:,:,i)=padarray(temp,[0,numberDiff/2,0],0);
        end
    end
    try
        diagMat2(:,:,i)=spdiags(fliplr(image));
    catch
        % Add a column at the beginning to prevent errors
        temp = spdiags(fliplr(image));
        numberDiff = abs(size(temp,2)-size(diagMat2,2));
        if mod(numberDiff,2) % Odd difference number
            temp = padarray(temp,[0,(numberDiff+1)/2,0],0);
            diagMat2(:,:,i) = temp(:,1:end-1);
        else
            diagMat2(:,:,i) = padarray(temp,[0,numberDiff/2,0],0);
        end
    end
end
for j = 1:szDiag(2)
    index = (diagMat1(:,j,1)~=0);
    nTemp = sum(index);
    image1 = zeros(sizeV(1),nTemp);
    image2 = zeros(sizeV(1),nTemp);
    for k = 1:sizeV(1)
        image1(k,1:nTemp) = diagMat1(index(1:end),j,k)';
        image2(k,1:nTemp) = diagMat2(index(1:end),j,k)';
    end
    uniqueIm = unique(image2);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image2;
    for i = 1:NLtemp
        indexRow(i) = find(uniqueIm(i)==levels);
        image2(temp==uniqueIm(i)) = i;
    end
    seq = zigzag(fliplr(image2));
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun);
end
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];
[textures] = getGLRLMtextures(GLRLM);
GLN = GLN + textures.GLN * sum(GLRLM(:)); % (Aerts et al., 2014) used the non-normalised version of GLN


GLN = GLN/13;
end