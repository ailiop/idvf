function [varargout] = mfigure (varargin)
%
% MFIGURE - Figure object with screen-adaptive size
%
% SYNTAX
%
%   MFIGURE 
%   MFIGURE( SCALE )
%   H = MFIGURE( ... )
%
% INPUT
%
%   SCALE       Figure-size to screen-size ratio        [scalar]
%               {default: 0.6}
%
% OUTPUT
%
%   H           Figure object handle                    [matlab.ui.figure]
%
% DESCRIPTION
%
%   MFIGURE is a wrapper around the FIGURE function.  It creates a figure
%   object whose size is set relative to the screen size.  This can be
%   useful when using large, high-resolution screens.
%
% DEPENDENCIES
%
%   <none>
%
%
% See also      figure
%
    
    
    % default value for optional scaling factor
    switch nargin
      case 0
        figureScale = 0.6;
      case 1
        figureScale = varargin{1};
    end

    % Get screen parameters
    graphicalRoot = groot; % Get groot handle
    screenSizeGroot = graphicalRoot.ScreenSize(3:4); % Get screen size
                                                     % width and height

    % Compute screen offset
    figSizeEdgeOffset = (1 - figureScale) * graphicalRoot.ScreenSize(3:4) / 2;

    % Open figure
    hf = figure; % Initialize figure

    % Get/set units
    figUnits = hf.Units; % Get current figure units (users may change defaults)
    hf.Units = graphicalRoot.Units; % Force units the same

    % Compute figure size in terms of width and height
    figSize = screenSizeGroot - figSizeEdgeOffset*2; % width, height offsets

    % Position and resize figure
    hf.Position = [figSizeEdgeOffset(1) figSizeEdgeOffset(2) ...
                   figSize(1) figSize(2)]; % left bottom width height

    % Set figure units back
    hf.Units = figUnits; 
    
    % optional output: figure handle
    if nargout > 0
        varargout{1} = hf;
    end
    
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Kevin Mattheus Moerman
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% RELEASE
%
%   1.0.2 - February 11, 2019
%
% CHANGELOG
%
%   1.0.2 (February 11, 2019) - Alexandros
%       + rudimentary documentation
%
%   0.0.1 (Dec 10, 2018) - Kevin
%       . initial version
%         [https://github.com/ailiop/idvf/issues/4#issue-389261860]
%
% ------------------------------------------------------------


