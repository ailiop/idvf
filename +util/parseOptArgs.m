function opt = parseOptArgs (dflt, varargin)
%
% PARSEOPTARGS - Optional arguments parsing
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
%   OPT = PARSEOPTARGS( DFLT, 'Name', Value, ... )
%
% INPUT
%
%   DFLT        Struct with default parameters          [struct]
%               (field names and values comprise the 
%               optional parameter name-value pairs)
%   <Input name-value pairs for non-default parameters> [varargin]
%
% OUTPUT
%
%   OPT         Input and default parameters            [struct]
%
% DESCRIPTION
%
%   OPT = PARSEOPTARGS(DFLT,'Name',Value,...) is a wrapper around the
%   MATLAB inputParser, to facilitate optional argument set-up and parsing.
%
%   If any input value is empty (e.g., []), then the default value is used.
%
% PARAMETERS
%
%   The following parameters are explicitly set for the underlying
%   inputParser.  Use a customized argument parser if these are not
%   appropriate.
%
%       CaseSensitive   = false
%       KeepUnmatched   = false
%       PartialMatching = true
%       StructExpand    = true
%
% DEPENDENCIES
%
%   <none>
%
%
% See also      inputParser
%
    
    
    %% INITIALIZATION
    
    ip = inputParser;
    
    ip.CaseSensitive   = false;
    ip.KeepUnmatched   = false;
    ip.PartialMatching = true;
    ip.StructExpand    = true;
    
    
    %% PARAMETERS
    
    argNames = fieldnames( dflt );
    for i = 1 : length(argNames)
        addParameter( ip, argNames{i}, dflt.(argNames{i}) );
    end
    
    
    %% PARSE AND RETURN
    
    parse( ip, varargin{:} );
    opt = ip.Results;
    
    
    %% SET EMPTY VALUES TO DEFAULTS
    
    for i = 1 : length(argNames)
        if isempty( opt.(argNames{i}) )
            opt.(argNames{i}) = dflt.(argNames{i});
        end
    end
    
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% RELEASE
%
%   1.0.0 - October 31, 2018
%
% ------------------------------------------------------------
