function [varargout] = checkFusionSliceDim(varargin)
% -------------------------------------------------------------------------
% function [varargout] = checkFusionSliceDim(varargin)
% -------------------------------------------------------------------------
% DESCRIPTION: 
% This function performs checks prior to PET/MRI fusion to ensure that the
% PET and MRI volumes have the same number of slices. Differences of one 
% slice sometimes occur (maybe due to interpolation in the MIM software).
%
% --> The sequence in the inputs and outputs needs to be PET, MRI, PET, MRI, etc.
% -------------------------------------------------------------------------
% REFERENCE:
% [1] Vallieres, M. et al. (2015). A radiomics model from joint FDG-PET and 
%     MRI texture features for the prediction of lung metastases in soft-tissue 
%     sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 
%     5471-5496. doi:10.1088/0031-9155/60/14/5471
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: May 2015
% -------------------------------------------------------------------------
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


if size(varargin{1},3) ~= size(varargin{2},3)
    if size(varargin{1},3) > size(varargin{2},3)
        varargout{2} = varargin{2};
        varargout{1} = varargin{1}(:,:,1:end-1);
        if nargin > 2
            for i = 3:nargin
                if mod(i,2)
                    varargout{i} = varargin{i}(:,:,1:end-1);
                else
                    varargout{i} = varargin{i};
                end
            end
        end
    else
        if abs(size(varargin{1},3) - size(varargin{2},3)) == 1
            varargout{1} = varargin{1};
            varargout{2} = varargin{2}(:,:,1:end-1);
            if nargin > 2
                for i = 3:nargin
                    if mod(i,2)
                        varargout{i} = varargin{i};
                    else
                        varargout{i} = varargin{i}(:,:,1:end-1);
                    end
                end
            end
        end
        if abs(size(varargin{1},3) - size(varargin{2},3)) == 2
            varargout{1} = varargin{1};
            varargout{2} = varargin{2}(:,:,2:end-1);
            if nargin > 2
                for i = 3:nargin
                    if mod(i,2)
                        varargout{i} = varargin{i};
                    else
                        varargout{i} = varargin{i}(:,:,2:end-1);
                    end
                end
            end
        end
    end
else
    for i = 1:nargin
        varargout{i} = varargin{i};
    end
end

end