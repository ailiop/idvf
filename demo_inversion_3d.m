% demo_inversion_3d.m
%
% 3D DVF inversion [1,2] demo script (3D-expanded 2D DVF).
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
% REFERENCES
%
%   [1] A. Dubey*, A.-S. Iliopoulos*, X. Sun, F.-F. Yin, and L. Ren,
%   "Iterative inversion of deformation vector fields with feedback
%   control," Medical Physics, vol. 45, no. 7, pp. 3147-3160, 2018.
%   DOI: 10.1002/mp.12962.
%
%   [2] A. Dubey, A.-S. Iliopoulos, X. Sun, F.-F. Yin, and L. Ren,
%   "Symmetric completion of deformation registration via bi-residual
%   inversion." In preparation.
%



%% ==================== CLEAN UP

clear variables
close all


%% ==================== PARAMETERS

% example DVF data
pathData   = 'data/c-deformation.mat';
ioData     = matfile( pathData );
F          = ioData.F;
[szDom, ~] = dvf.sizeVf( F );

% 2D-to-3D DVF embedding (illustration with zero Z-axis displacement)
zdim = 15;
F = repmat( F, [1 1 1 zdim] );          % (NX x NY x 2 x NZ)
F = permute( F, [1 2 4 3] );            % (NX x NY x NZ x 2)
F = cat( 4, F, zeros( [szDom zdim] ) ); % (NX x NY x NZ x 3)

% inversion parameters
opt.control       = 'midrange';
opt.scale         = 1;
opt.numIteration  = 20;
opt.stopThreshold = 1e-6;
opt.accThreshold  = 0;
opt.InitialValue  = [];
opt.Mask          = dvf.maskDomain( F );

% create an array of inversion parameters to compare multiple results
opt = repmat( opt, [7 1] );
opt(1).control       = 0.0;             % constant mu = 0
opt(2).control       = 0.5;             % constant mu = 0.5
opt(3).control       = 'alternating';   % alternating
opt(4).control       = 'optimal';       % locally optimal
opt(5).control       = 'midrange';      % local midrange
opt(6).control       = 'midrange';      % local midrange & acc.
opt(6).accThreshold  = 1;               %  .
opt(7).control       = 'midrange';      % local midrange & acc. & multiscale
opt(7).accThreshold  = 1;               %  .
opt(7).scale         = [0.5, 1.0];      %  .
opt(7).numIteration  = [4 16];          %  .

% legend labels for each control setting
lgnd = {'\mu = 0', '\mu = 0.5', '\mu_{oe}', '\mu_*(x)', '\mu_m(x)', ...
        '\mu_m(x) & acc', '\mu_m(x) & acc & MS'};

% IC residual percentile curve visualization options
prcIcMeasure  = [50 90 95 98 100];
argVisIcCurve = {'-o', 'LineWidth', 2};


%% ==================== (BEGIN)

fprintf( '\n***** BEGIN (%s) *****\n\n', mfilename );


%% ==================== VISUALIZATION: REFERENCE & DEFORMED IMAGES

fprintf( '...visualizing synthetic reference and deformed images...\n' );

IRef = util.imageGrid( [szDom zdim] );
IDfm = dvf.imdeform( IRef, F );

figure

% reference image
subplot( 1, 2, 1 )
imagesc( IRef(:,:,(zdim+1)/2).' )
axis image
colormap( gray )
title( {'reference image', '(central axial slice)'} )
drawnow

% deformed (study) image
subplot( 1, 2, 2 )
imagesc( IDfm(:,:,(zdim+1)/2).' )
axis image
colormap( gray )
title( {'deformed image', '(central axial slice)'} )
drawnow


%% ==================== VISUALIZATION: SPECTRAL NTDC MEASURES

fprintf( '...calculating & visualizing spectral NTDC measures...\n' );

[CtrlIdx, Rho, Det, Lambda] = dvf.ntdcMeasures( F );

figure

% control index
subplot( 1, 3, 1 )
imagesc( CtrlIdx(:,:,(zdim+1)/2).' )
axis image
colormap( jet )
colorbar
title( {'algebraic control index', '(central axial slice)'} );
drawnow

% NTDC spectral radius
subplot( 1, 3, 2 )
imagesc( Rho(:,:,(zdim+1)/2).' )
axis image
colormap( jet )
colorbar
title( {'NTDC spectral radius', '(central axial slice)'} );
drawnow

% determinant map
subplot( 1, 3, 3 )
imagesc( Det(:,:,(zdim+1)/2).' )
axis image
colormap( jet )
colorbar
title( {'determinant map', '(central axial slice)'} );
drawnow


%% ==================== INVERSION

fprintf( '...computing inverse with %d iteration schemes...\n', length(opt) );

G      = cell( size(opt) );
prctRG = cell( size(opt) );
prctRF = cell( size(opt) );

for i = 1 : length(opt)
    fprintf( '   - scheme #%d (%s)...\n', i, lgnd{i} );
    opt(i).control = dvf.feedbackControlVal( Lambda, opt(i).control );
    [G{i}, ~, ~, prctRG{i}, prctRF{i}] = dvf.inversion( F, opt(i) );
end


%% ==================== VISUALIZATION: IC RESIDUALS PER ITERATION

fprintf( '...visualizing IC residual & error percentile curves...\n' );

% resolution-normalized iteration steps
kstep = cell( size(opt) );
for i = 1 : length(opt)
    kstep{i} = arrayfun( @(s,k) repmat( s^2 + (s~=1)/k, [1 k] ), ...
                         opt(i).scale, opt(i).numIteration, ...
                         'UniformOutput', false );
    kstep{i} = cumsum( horzcat( 0, kstep{i}{:} ) );
end

figure
numPrct = length( prcIcMeasure );

% legend
subplot( 5, numPrct, (1:numPrct) )
plot( nan( 2, length(opt) ), argVisIcCurve{:} );
title( 'step-wise IC residual percentile curves' )
legend( lgnd, 'Orientation', 'horizontal', 'Location', 'south' )
axis( gca, 'off' )
drawnow

for p = 1 : numPrct                     % ---- each percentile
    
    % IC residual (RG)
    subplot( 5, numPrct, p + (1:2)*numPrct )
    hold on
    for i = 1 : length(opt)
        plot( kstep{i}, prctRG{i}(p,:), argVisIcCurve{:} );
    end
    ylabel( sprintf( 'r_G[%d%%]', prcIcMeasure(p) ) );
    xlabel( 'amortizerd iteration step' )
    set( gca, 'YScale', 'log' )
    grid on
    drawnow
    
    % IC residual (RF)
    subplot( 5, numPrct, p + (3:4)*numPrct )
    hold on
    for i = 1 : length(opt)
        plot( kstep{i}, prctRF{i}(p,:), argVisIcCurve{:} );
    end
    ylabel( sprintf( 'r_F[%d%%]', prcIcMeasure(p) ) );
    xlabel( 'amortized iteration step' )
    set( gca, 'YScale', 'log' )
    grid on
    drawnow
    
end  % for (p)


%% ==================== (END)

fprintf( '\n***** END (%s) *****\n\n', mfilename );



%%------------------------------------------------------------
%
% AUTHORS
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% VERSION
%
%   0.1 - October 08, 2018
%
% CHANGELOG
%
%   0.1 (Oct 08, 2018) - Alexandros
%       * initial implementation
%
% ------------------------------------------------------------
