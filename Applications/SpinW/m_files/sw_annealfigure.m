function hFigure = sw_annealfigure
% creates a figure for displaying the status of the annealing simulation
%
% hFigure = SW_ANNEALFIGURE()
%
% See also SW.ANNEAL.
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

% Position of the new figure window.
posFig = get(0,'DefaultFigurePosition');
posFig = [posFig(1:2) 400 400];

% Create new figure.
hFigure = figure;
set(0,'Showhidden','on')

set(hFigure,...
    'Position',      posFig,...
    'DockControls',  'off',...
    'PaperType',     'A4',...
    'Tag',           'sw_anneal',...
    'Toolbar',       'figure',...
    'Name',          'SpinW : Annealing status');

end