function [Vnew] = imresize3D(V,scale,tsize,ntype,npad)
% -------------------------------------------------------------------------
% function [Vnew] = imresize3D(V,scale,tsize,ntype,npad)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function resizes a 3D input volume to new dimensions. It is adapted 
% from the original code of D. Kroon of Twente(July 2008) available at: 
% <http://www.mathworks.com/matlabcentral/fileexchange/21451-multimodality-non-rigid-demon-algorithm-image-registration/content//functions/imresize3d.m>
% -------------------------------------------------------------------------
% INPUTS:
% - V: The input volume (3D array).
% - scale: scaling factor, when used set tsize to [].
% - nsize: new dimensions, when used set scale to [].
% - ntype: Type of interpolation ('nearest', 'linear', or 'cubic')
% - npad: Boundary condition ('replicate', 'symmetric', 'circular', 'fill', 
%         or 'bound')  
% -------------------------------------------------------------------------
% OUTPUTS:
% - Vnew: Resized volume.
% -------------------------------------------------------------------------
% AUTHOR(S): 
% - Dirk-Jan Kroon
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: July 2008 (D. Kroon)
% - Revision: May 2015 (M. Vallieres)
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
%
%    _______________________________________________________________
%
% --> Copyright (c) 2009, Dirk-Jan Kroon
%     All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% -------------------------------------------------------------------------


% Check the inputs
if(exist('ntype', 'var') == 0), ntype='nearest'; end
if(exist('npad', 'var') == 0), npad='bound'; end
if(exist('scale', 'var')&&~isempty(scale)), tsize=round(size(V)*scale); end
if(exist('tsize', 'var')&&~isempty(tsize)),  scale=(tsize./size(V)); end

% Make transformation structure   
T = makehgtform('scale',scale);
tform = maketform('affine', T);

% Specify resampler
R = makeresampler(ntype, npad);

% Resize the image volume
Vnew = tformarray(V, tform, R, [1 2 3], [1 2 3], tsize, [], 0);
end