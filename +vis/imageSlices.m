function imageSlices (I, varargin)
%
% IMAGESLICES - Slice visualization of 3D image in subplots
%
% ----------------------------------------------------------------------
%
% Copyright (C) 2018, Department of Computer Science, Duke University
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% ----------------------------------------------------------------------
%   
% SYNTAX
%
%   IMAGESLICES( I )
%   IMAGESLICES( I, 'Name', Value, ... )
%
% INPUT
%
%   I           3D image                        [NX x NY x NZ]
%
% OPTIONS
%
%   'Slices'            [3-vector]              {ceil( [NX NY NZ] / 2 )}
%
%       Image slice indices (along X, Y, and Z axes).
%
%   'Colormap'          [C-by-3]                {gray}
%
%       IMAGESC colormap.
%
%   'CAxis'             [2-vector | 'auto']     {'auto'}
%
%       Color-axis limits (input to CAXIS) function.
%
%   'Rows'              [scalar]                {1}
%
%       Number of subplot rows.
%
%   'Columns'           [scalar]                {3}
%
%       Number of subplot columns.
%
%   'StartIndex'        [scalar]                {1}
%
%       First subplot index.  The slices are visualized in Rows x Columns
%       subplots, from index StartIndex (XY/axial slice), to indices
%       StartIndex+1 (XZ/coronal slice) and StartIndex+2 (YZ/sagittal slice).
%
%   'Title'             [string]                {''}
%
%       Title string to be shared among slice subplots.
%
%   'PColor'            [boolean]               {false}
%
%       If TRUE, then IMAGESLICES visualizes the image slices using PCOLOR
%       instead of IMAGESC.  This places a little more strain on the
%       renderer, but can be useful for distinguishing NaN values.
%
% DESCRIPTION
%
%   IMAGESLICES is a convenience function for visualizing 3 dimension-wise
%   slices of a 3D image in subplots.
%
% DEPENDENCIES
%
%   util.parseOptArgs
%
%
% See also      imagesc, pcolor, slice, volumeViewer
%
    
    
    %% PARAMETERS
    
    % options and default values
    opt.slices     = ceil( size(I) / 2 );
    opt.colormap   = gray;
    opt.caxis      = 'auto';
    opt.rows       = 1;
    opt.columns    = 3;
    opt.startIndex = 1;
    opt.title      = '';
    opt.pcolor     = false;
    
    % parse optional arguments
    opt = util.parseOptArgs( opt, varargin{:} );
    
    
    %% IMAGESC vs PCOLOR
    
    if opt.pcolor
        hPlot = @pcolor;
    else
        hPlot = @imagesc;
    end
    
    
    %% VISUALIZATION
    
    % XY slice
    subplot( opt.rows, opt.columns, opt.startIndex )
    hPlot( squeeze( I(:,:,opt.slices(3)) ).' )
    axis image
    shading flat
    caxis( opt.caxis )
    colormap( opt.colormap )
    colorbar
    title( {opt.title, sprintf( 'XY slice (z = %d)', opt.slices(3) )} )
    drawnow
    
    % XZ slice
    subplot( opt.rows, opt.columns, opt.startIndex+1 )
    hPlot( squeeze( I(:,opt.slices(2),:) ).' )
    axis image
    shading flat
    caxis( opt.caxis )
    colormap( opt.colormap )
    colorbar
    title( {opt.title, sprintf( 'XZ slice (y = %d)', opt.slices(2) )} )
    drawnow
    
    % YZ slice
    subplot( opt.rows, opt.columns, opt.startIndex+2 )
    hPlot( squeeze( I(opt.slices(1),:,:) ).' )
    axis image
    shading flat
    caxis( opt.caxis )
    colormap( opt.colormap )
    colorbar
    title( {opt.title, sprintf( 'YZ slice (x = %d)', opt.slices(1) )} )
    drawnow
    
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% RELEASE
%
%   1.0.2 - February 11, 2019
%
% CHANGELOG
%
%   1.0.2 (Feb 11, 2019) - Alexandros
%       . initial version
%
% ------------------------------------------------------------
