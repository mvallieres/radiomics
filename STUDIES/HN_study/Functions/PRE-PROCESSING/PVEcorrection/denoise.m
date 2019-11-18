function volumeOut = denoise(volumeIn,level)
% -------------------------------------------------------------------------
% function volumeOut = denoise(volumeIn,level)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function uses a wavelet-based thresholding approach to denoise an 
% input volume. The premise of the program is to apply the wavelet 
% transform to a volume, apply the thresholding operator to the transformed 
% volume, and finally apply an inverse wavelet transform to obtain the 
% denoised volume. As we need translation invariance, we use the
% nondecimated wavelet transform (ndwt). This is limited to 2D, thus we 
% operate on all three planes indiviudallly, followed by averaging on a 
% per-voxel basis. See ref. [1] for more details.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Boussion, N. et al. (2009). Incorporation of wavelet-based denoising
%     in iterative deconvolution for partial volume correction in 
%     whole-body PET imaging. Eur J Nucl Med Mol Imaging, 36(7), 1064-1075.
% -------------------------------------------------------------------------
% INPUTS:
% - volumeIn: 3D array representing the input PET volume to correct for PVE.
% - level: Numeric value specifying the level of decmposition of the
%          wavelet transform. It is recommended to use level = 3.
% -------------------------------------------------------------------------
% OUTPUTS:
% - volumeOut: 3D array representing the denoised input volume.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Andre Diamant <adboustead@gmail.com>
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2016
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Martin Vallieres
%
%    This package is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This package is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this package.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

global wavelet_name

% FIRST PLANE
counter=1;
index_plane_1=cell([1,size(volumeIn,3)]);
denoisedplane_1=cell([1,size(volumeIn,3)]);   %Initiate all our cells/arrays
resultvector_1=zeros(size(volumeIn));

for i=1:size(volumeIn,3), %The "plane" is in 2-d and the "i" represents the 3rd dimension. Each "i" component is a different slice of the plane.
    
    results=cell([3*level+1,1]); %Create a "results" cell to place the denoised results
    plane_1=volumeIn(:,:,i);
    index_plane_1{1,i}=zeros(size(volumeIn,1),size(volumeIn,2));            %Set our index plane with indices in increasing order.
    index_plane_1{1,i}(1:end)=counter:(counter-1+numel(index_plane_1{1,i}));
    counter=counter+numel(index_plane_1{1,i});      %Counter required to ensure that subsequent slices start at the correct index
    
    for j=1:level,
        %Then we apply the wavelet transform to our plane (must be 2-d) (up
        %to scale j)

        if ~isempty(wavelet_name)
            wavelet = ndwt2(plane_1,j,wavelet_name);
        else
            wavelet = ndwt2(plane_1,j,'bior3.5'); 
        end 

        %Second we apply the threshold operator T to the wavelet space

        columnwavelet = wavelet.dec{1}(:); % Turn first subband into a column

        M=median(abs(columnwavelet)); %Median of the first subband
        noisevar=(M/0.6745)^2;
        var=sum(columnwavelet.^2)/(size(columnwavelet,2).^2);  %Var and noisevar are chosen by the BayesShrink method

        if var-noisevar <= 0,
            threshold = 1000000000; %This could be infinity, simply to discount daa if var < noisevar. (Shouldn't happen)
            disp('Threshold set to big')

        else    
            threshold=noisevar./sqrt(var-noisevar);
        end  


        if j==level
            for x=1:4,
            results{x}=sign(wavelet.dec{x}).*max(0, abs(wavelet.dec{x})-threshold);
            end
            
        else
            for x=2:4,
            results{3*(level-j)+x}=sign(wavelet.dec{x}).*max(0, abs(wavelet.dec{x})-threshold); %Actual operator T applied, with a threshold determined as above
            end
        end
                     
    end

    wavelet.dec=results;
    
   %Finally, apply inverse wavelet transform to get it back
      
    denoisedplane_1{i}=indwt2(wavelet);
    
    resultvector_1(index_plane_1{1,i}(:))=denoisedplane_1{i};
    
end

%disp('Plane 1 done!')


% SECOND PLANE
index_plane_2=cell([1,size(volumeIn,2)]);
denoisedplane_2=cell([1,size(volumeIn,2)]);       %Initiate all our cells/arrays
resultvector_2=zeros(size(volumeIn));


for i=1:size(volumeIn,2),     %The "plane" is in 2-d and the "i" represents the 3rd dimension. Each "i" component is a different slice of the plane.
    
    results=cell([3*level+1,1]); %Create a "results" cell to place the denoised results
    plane_2=zeros(size(volumeIn,1),size(volumeIn,3));
    plane_2(:,1:end)=volumeIn(:,i,1:end);
    index_plane_2{1,i}=zeros(size(volumeIn,1),size(volumeIn,3));
    
    for j=1:size(volumeIn,3),
        index_plane_2{1,i}(:,j)=index_plane_1{1,j}(:,i);
    end
        
    for j=1:level,
        %Then we apply the wavelet transform to our inputed area (must be 2-d)
        
        if ~isempty(wavelet_name)
            wavelet = ndwt2(plane_2,j,wavelet_name);
        else
            wavelet = ndwt2(plane_2,j,'bior3.5'); 
        end
        
        %Second we apply the threshold operator T to the wavelet space

        columnwavelet = wavelet.dec{1}(:); % Turn first subband into a column

        M=median(abs(columnwavelet)); %Median of the first subband
        noisevar=(M./0.6745)^2;
        var=sum(columnwavelet.^2)/(size(columnwavelet,2).^2);

        if var-noisevar <= 0,
            threshold = 1000000000; %This could be infinity
            disp('Threshold set to big')

        else    
            threshold=noisevar./sqrt(var-noisevar);
%            fprintf('Threshold good! = %8.8f \n', threshold)
        end  


        if j==level
            for x=1:4,
            results{x}=sign(wavelet.dec{x}).*max(0, abs(wavelet.dec{x})-threshold); %Actual operator T applied
            end
            
        else
            for x=2:4,
            results{3*(level-j)+x}=sign(wavelet.dec{x}).*max(0, abs(wavelet.dec{x})-threshold); %Actual operator T applied
            end
        end
                     
    end

    wavelet.dec=results;
    
   %Finally, apply inverse wavelet transform to get it back
      
    denoisedplane_2{i}=indwt2(wavelet);
   
    resultvector_2(index_plane_2{1,i}(:))=denoisedplane_2{i};
    
end

%disp('Plane 2 done!')

% THIRD PLANE

index_plane_3=cell([1,size(volumeIn,1)]);
denoisedplane_3=cell([1,size(volumeIn,1)]);
resultvector_3=zeros(size(volumeIn));

for i=1:size(volumeIn,1),
    
    results=cell([3*level+1,1]); %Create a "results" cell to place the denoised results
    plane_3=zeros(size(volumeIn,2),size(volumeIn,3));
    plane_3(:,1:end)=volumeIn(i,:,1:end);
    index_plane_3{1,i}=zeros(size(volumeIn,2),size(volumeIn,3));
    
    for j=1:size(volumeIn,3),
        index_plane_3{1,i}(:,j)=index_plane_1{1,j}(i,:);
    end
    
    for j=1:level,
        %Then we apply the wavelet transform to our inputed area (must be 2-d)

        if ~isempty(wavelet_name)
            wavelet = ndwt2(plane_3,j,wavelet_name);
        else
            wavelet = ndwt2(plane_3,j,'bior3.5'); 
        end 
            
        %Second we apply the threshold operator` T to the wavelet space

        columnwavelet = wavelet.dec{1}(:); % Turn first subband into a column

        M=median(abs(columnwavelet)); %Median of the first subband
        noisevar=(M./0.6745)^2;
        var=sum(columnwavelet.^2)/(size(columnwavelet,2).^2);

        if var-noisevar <= 0,
            threshold = 1000000000; %This could be infinity
            disp('Threshold set to big')

        else    
            threshold=noisevar./sqrt(var-noisevar);
%             fprintf('Threshold good! = %8.8f \n', threshold)
        end  


        if j==level
            for x=1:4,
            results{x}=sign(wavelet.dec{x}).*max(0, abs(wavelet.dec{x})-threshold); %Actual operator T applied
            end
            
        else
            for x=2:4,
            results{3*(level-j)+x}=sign(wavelet.dec{x}).*max(0, abs(wavelet.dec{x})-threshold); %Actual operator T applied
            end
        end
                     
    end

    wavelet.dec=results;
    
   %Finally, apply inverse wavelet transform to get it back
      
    denoisedplane_3{i}=indwt2(wavelet);

    resultvector_3(index_plane_3{1,i}(:))=denoisedplane_3{i};

    
end

%disp('Plane 3 done!')

denoisedvolumevector=(resultvector_1+resultvector_2+resultvector_3)/3;      %Take the mean of our three results

denoisedvolume=reshape(denoisedvolumevector,[size(volumeIn,1),size(volumeIn,2),size(volumeIn,3)]);

volumeOut = denoisedvolume;         %Return our denoised volume!

end

    