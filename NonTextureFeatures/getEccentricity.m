function [eccentricity] = getEccentricity(ROIonly,pixelW,sliceS)
% -------------------------------------------------------------------------
% function [eccentricity] = getEccentricity(ROIonly,pixelW,sliceS)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function computes the eccentricity metric of the region of interest 
% (ROI) of an input volume. The ellipsoid that best fits the ROI is first 
% computed following the framework of reference [1], using an adaptation of
% a code available at: <http://www.mathworks.com/matlabcentral/fileexchange/23377-ellipsoid-fitting>
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Li, Q. and Griffiths, J. G. (2004). Least-square ellipsoid specific
%     fitting. Proceedings of the Geometric Modeling and Processing (GMP), 
%     Beijing, China, 335-340.
% -------------------------------------------------------------------------
% INPUTS:
% - ROIonly: 3D array, with voxels outside the ROI set to NaNs.
% - pixelW: Pixel width, or in-plane resolution, in mm.
% - sliceS: Slice spacing, in mm.
% -------------------------------------------------------------------------
% OUTPUTS:
% - eccentricity: Metric given by sqrt[1 - a*b/c^2] where c is the longest 
%                 semi-principal axes of the ellipsoid, and a and b are the 
%                 second and third longest semi-principal axes of the 
%                 ellipsoid.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Q. Li
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: March 2009 (Q. Li)
% - Adaptation: May 2015 (M. Vallieres)
%--------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Q. Li, Martin Vallieres
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

mask = ~isnan(ROIonly); % Find mask covering the ROI

% ISOTROPIC RESAMPLING
sFactor = sliceS/pixelW; % scaling factor
mask = imresize3D(mask,[],[round(double(size(mask,1))),round(double(size(mask,2))),round(double(size(mask,3))*sFactor)],'nearest','fill');


% INITIALIZATION OF ECCENTRICITY COMPUTATION
perimeter = bwperim(mask,18);
nPoints = length(find(perimeter));
X = zeros(nPoints,1); Y = zeros(nPoints,1); Z = zeros(nPoints,1);
count = 1;
for i = 1:size(perimeter,3)
    [row,col] = find(perimeter(:,:,i));
    p = length(row);
    if p > 0
        X(count:count+p-1,1) = col(1:end);
        Y(count:count+p-1,1) = row(1:end);
        Z(count:count+p-1,1) = i;
        count = count + p;
    end
end


% LI'S ELLIPSOID FITTING ALGORITHM
dx = X(:); dy = Y(:); dz = Z(:);
n = size(dx,1);
D = [dx.*dx, dy.*dy,  dz.*dz, 2.*dy.*dz, 2.*dx.*dz, 2.*dx.*dy, ...
        2.*dx, 2.*dy, 2.*dz, ones(n,1)]';
S = D*D';

% Create constraint matrix C:
C(6,6)=0; C(1,1)=-1; C(2,2)=-1; C(3,3)=-1; C(4,4)=-4; C(5,5)=-4; C(6,6)=-4;
C(1,2)=1; C(2,1)=1; C(1,3)=1; C(3,1)=1; C(2,3)=1; C(3,2)=1;

% Solve generalized eigensystem
S11 = S(1:6, 1:6); S12 = S(1:6, 7:10); S22 = S(7:10,7:10);
A = S11-S12*pinv(S22)*S12'; CA = inv(C)*A;
[gevec, geval] = eig(CA);

% Find the largest eigenvalue(the only positive eigenvalue)
In = 1;
maxVal = geval(1,1);
for i = 2:6
   if (geval(i,i) > maxVal)
       maxVal = geval(i,i);
       In = i;
   end;
end;

% Find the fitting
v1 = gevec(:, In); 
v2 = -pinv(S22)*S12'*v1;
v = [v1; v2];

% Algebraic from of the ellipsoid
A = [ v(1) v(6) v(5) v(7); ...
  v(6) v(2) v(4) v(8); ...
  v(5) v(4) v(3) v(9); ...
  v(7) v(8) v(9) -1 ];

% Center of the ellipsoid
center = -A( 1:3, 1:3 ) \ [ v(7); v(8); v(9) ];

% Corresponding translation matrix
T = eye( 4 );
T(4, 1:3 ) = center';

% Translating to center
R = T * A * T';

% Solving eigenproblem
[ evecs evals ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
radii = sqrt( 1 ./ diag( evals ) );


% ECCENTRICITY COMPUTATION
eccentricity = sqrt(1 - (radii(2)*radii(3)/radii(1)^2));

end