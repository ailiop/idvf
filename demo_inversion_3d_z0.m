% demo_inversion_3d_z0.m
%
% 3D DVF inversion [1,2] demo script (3D-expanded 2D DVF; no displacement
% along Z-axis).
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
%   - DOI: 10.1002/mp.12962
%   - arXiv: 1610.08589 [cs.CV]
%
%   [2] A. Dubey, "Symmetric completion of deformable registration via
%   bi-residual inversion," PhD thesis, Duke University, Durham, NC, USA,
%   2018.
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
zdim  = 15;
F     = repmat( F, [1 1 1 zdim] );          % (NX x NY x 2 x NZ)
F     = permute( F, [1 2 4 3] );            % (NX x NY x NZ x 2)
F     = cat( 4, F, zeros( [szDom zdim] ) ); % (NX x NY x NZ x 3)
szDom = [szDom zdim];

% inversion parameters
opt.control       = 'midrange';
opt.scale         = 1;
opt.numIteration  = 20;
opt.stopThreshold = 1e-6;
opt.accThreshold  = 0;
opt.InitialValue  = [];
opt.Mask          = dvf.maskDomain( F );

% create an array of inversion parameters to compare multiple schemes
opt = repmat( opt, [8 1] );
opt(1).control       = 0.0;             % constant mu = 0
opt(2).control       = 0.5;             % constant mu = 0.5
opt(3).control       = 'alternating';   % alternating
opt(4).control       = 'optimal';       % pointwise optimal
opt(5).control       = 'midrange';      % local midrange
opt(6).control       = 'midrange';      % local midrange & acc.
opt(6).accThreshold  = 1;               %  .
opt(7).control       = 'midrange';      % local midrange & acc. & multiscale
opt(7).accThreshold  = 1;               %  .
opt(7).scale         = [0.5, 1.0];      %  .
opt(7).numIteration  = [4 16];          %  .
opt(8).control       = 'optimal';       % pointwise optimal & acc. & multiscale
opt(8).accThreshold  = 1;               %  .
opt(8).scale         = [0.5, 1.0];      %  .
opt(8).numIteration  = [4 16];          %  .

% legend labels for each control setting
lgnd = {'\mu = 0', '\mu = 0.5', '\mu_{oe}', '\mu_*(x)', '\mu_m(x)', ...
        '\mu_m(x) & acc', '\mu_m(x) & acc & MS', '\mu_*(x) & acc & MS'};

% IC residual percentile curve visualization options
prcIcMeasure  = [50 90 95 98 100];
argVisIcCurve = {'-o', 'LineWidth', 2};
clrIcCurve    = [0.8650 0.8110 0.4330;  % yellow
                 0.9718 0.5553 0.7741;  % pink
                 0.6859 0.4035 0.2412;  % brown
                 1.0000 0.5482 0.1000;  % orange
                 0.6365 0.3753 0.6753;  % purple
                 0.3718 0.7176 0.3612;  % green
                 0.2941 0.5447 0.7494;  % blue
                 0.9047 0.1918 0.1988]; % red

% IC residual maps visualization options
zIdx      = (zdim + 1) / 2;             % axial slice index
climIcMgn = [0 1];                      % color-axis limits


%% ==================== (BEGIN)

fprintf( '\n***** BEGIN (%s) *****\n\n', mfilename );


%% ==================== VISUALIZATION: REFERENCE & STUDY IMAGES

fprintf( '...visualizing synthetic reference and study images...\n' );

IRef = util.imageGridSmooth( szDom, [15 15 5] );
IStd = dvf.imdeform( IRef, F );

hFig      = vis.mfigure;
hFig.Name = 'reference and study (deformed) images';

% reference image
subplot( 1, 2, 1 )
imagesc( IRef(:,:,zIdx).' )
axis image
colormap( gray )
title( {'reference image', sprintf( '(axial slice z = %d)', zIdx )} )
drawnow

% study (deformed) image
subplot( 1, 2, 2 )
imagesc( IStd(:,:,zIdx).' )
axis image
colormap( gray )
title( {'study image', sprintf( '(axial slice z = %d)', zIdx )} )
drawnow


%% ==================== VISUALIZATION: SPECTRAL NTDC MEASURES

fprintf( '...calculating & visualizing spectral NTDC measures...\n' );

[CtrlIdx, Rho, Det, Lambda] = dvf.ntdcMeasures( F );

hFig      = vis.mfigure;
hFig.Name = 'spectral NTDC measures of forward DVF';

% control index
subplot( 1, 3, 1 )
imagesc( CtrlIdx(:,:,zIdx).' )
axis image
colormap( jet )
colorbar
title( {'algebraic control index', sprintf( '(axial slice z = %d)', zIdx )} );
drawnow

% NTDC spectral radius
subplot( 1, 3, 2 )
imagesc( Rho(:,:,zIdx).' )
axis image
colormap( jet )
colorbar
title( {'NTDC spectral radius', sprintf( '(axial slice z = %d)', zIdx )} );
drawnow

% determinant map
subplot( 1, 3, 3 )
imagesc( Det(:,:,zIdx).' )
axis image
colormap( jet )
colorbar
title( {'determinant map', sprintf( '(axial slice z = %d)', zIdx )} );
drawnow


%% ==================== INVERSION

fprintf( '...computing inverse with %d iteration schemes...\n', length(opt) );

G      = cell( size(opt) );
prctRG = cell( size(opt) );
prctRF = cell( size(opt) );
RG     = cell( size(opt) );
RF     = cell( size(opt) );

for i = 1 : length(opt)
    fprintf( '   - scheme #%d (%s)...\n', i, lgnd{i} );
    opt(i).control = dvf.feedbackControlVal( Lambda, opt(i).control );
    [G{i}, ~, ~, prctRG{i}, prctRF{i}] = dvf.inversion( F, opt(i) );
    RG{i} = dvf.icResidual( G{i}, F );
    RF{i} = dvf.icResidual( F, G{i} );
end


%% ==================== VISUALIZATION: IC RESIDUALS PER ITERATION

fprintf( '...visualizing step-wise IC residual percentile curves...\n' );

% resolution-normalized iteration steps
kstep = cell( size(opt) );
for i = 1 : length(opt)
    kstep{i} = arrayfun( @(s,k) repmat( s^2 + (s~=1)/k, [1 k] ), ...
                         opt(i).scale, opt(i).numIteration, ...
                         'UniformOutput', false );
    kstep{i} = cumsum( horzcat( 0, kstep{i}{:} ) );
end

hFig      = vis.mfigure;
hFig.Name = 'step-wise IC residual percentile curves';
numPrct   = length( prcIcMeasure );

% legend
subplot( 5, numPrct, (1:numPrct), ...
         'NextPlot', 'ReplaceChildren', 'ColorOrder', clrIcCurve )
plot( nan( 2, length(opt) ), argVisIcCurve{:} );
title( 'step-wise IC residual percentile curves' )
legend( lgnd, 'Orientation', 'horizontal', 'Location', 'south' )
axis( gca, 'off' )
drawnow

for p = 1 : numPrct                     % ---- each percentile
    
    % study IC residuals
    subplot( 5, numPrct, p + (1:2)*numPrct )
    hold on
    for i = 1 : length(opt)
        plot( kstep{i}, prctRG{i}(p,:), argVisIcCurve{:}, ...
              'Color', clrIcCurve(i,:) );
    end
    ylabel( sprintf( 's_k[%d%%]', prcIcMeasure(p) ) );
    xlabel( 'amortizerd iteration step (k)' )
    set( gca, 'YScale', 'log' )
    grid on
    drawnow
    
    % reference IC residuals
    subplot( 5, numPrct, p + (3:4)*numPrct )
    hold on
    for i = 1 : length(opt)
        plot( kstep{i}, prctRF{i}(p,:), argVisIcCurve{:}, ...
              'Color', clrIcCurve(i,:) );
    end
    ylabel( sprintf( 'r_k[%d%%]', prcIcMeasure(p) ) );
    xlabel( 'amortized iteration step (k)' )
    set( gca, 'YScale', 'log' )
    grid on
    drawnow
    
end  % for (p)


%% ==================== VISUALIZATION: IC RESIDUAL MAPS

fprintf( '...visualizing final IC residual maps...\n' );
fprintf( '   - white pixels indicate NaN (out-of-bounds) values\n' );

% calculate point-wise IC residual magnitudes
RGMgn = cellfun( @(r) sqrt( sum( r.^2, 4 ) ), RG, 'UniformOutput', false );
RFMgn = cellfun( @(r) sqrt( sum( r.^2, 4 ) ), RF, 'UniformOutput', false );

% study IC residuals
hFig      = vis.mfigure;
hFig.Name = 'study IC residual maps';
for i = 1 : length(opt)
    subplot( 2, 4, i )
    pcolor( RGMgn{i}(:,:,zIdx).' )
    axis image
    shading flat
    colormap jet
    caxis( climIcMgn )
    colorbar
    title( {sprintf( 'study IC residuals (%s)', lgnd{i} ), ...
            sprintf( '(axial slice z = %d)', zIdx )} )
    drawnow
end

% reference IC residuals
hFig      = vis.mfigure;
hFig.Name = 'reference IC residual maps';
for i = 1 : length(opt)
    subplot( 2, 4, i )
    pcolor( RFMgn{i}(:,:,zIdx).' )
    axis image
    shading flat
    colormap jet
    caxis( climIcMgn )
    colorbar
    title( {sprintf( 'reference IC residuals (%s)', lgnd{i} ), ...
            sprintf( '(axial slice z = %d)', zIdx )} )
    drawnow
end


%% ==================== VISUALIZATION: RECOVERED REFERENCE IMAGE

fprintf( '...reference image recovery...\n' );
fprintf( '   - I_rec(x) := I_std(x + v(x))\t' );
fprintf( ' [I_std(y) = I_ref(y + u(y))]\n' );
fprintf( '               = I_ref(x + v(x) + u(x + v(x)))\n' );
fprintf( '               = I_ref(x + s(x))\n' );

% ---------- calculation
fprintf( '   - calculation...\n' );

IRefRec = cell( size(opt) );
for i = 1 : length(opt)
    IRefRec{i} = dvf.imdeform( IRef, RG{i} );
end

% ---------- visualization
fprintf( '   - difference image visualization (I_ref - I_rec)...\n' );

hFig      = vis.mfigure;
hFig.Name = 'image-space errors in reference image recovery';
for i = 1 : length(opt)
    subplot( 2, 4, i )
    pcolor( (IRef(:,:,zIdx) - IRefRec{i}(:,:,zIdx)).' );
    axis image
    shading flat
    caxis( [-1, +1] )
    colormap parula
    colorbar
    title( {lgnd{i}, 'I_{ref}(x) - I_{ref}(x+s(x))'} )
    drawnow
end


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
%   0.2 - December 21, 2018
%
% CHANGELOG
%
%   0.2 (Dec 21, 2018) - Alexandros
%       + added IC residual magnitude maps
%       + added image-space error maps (reference image recovery)
%       + added control scheme: pointwise optimal control values with
%         local search, local acceleration, and two-scale iteration
%       . changed synthetic reference image to smooth version
%       . parameter clean-up and explicit visualization options
%
%   0.1 (Oct 08, 2018) - Alexandros
%       . initial implementation
%
% ------------------------------------------------------------
