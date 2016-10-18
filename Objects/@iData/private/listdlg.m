function [selection,value] = listdlg_nonmodal(varargin)
%LISTDLG_nonmodal  List selection dialog box.
%   [SELECTION,OK] = LISTDLG('ListString',S) creates a modal dialog box
%   which allows you to select a string or multiple strings from a list.
%   SELECTION is a vector of indices of the selected strings (length 1 in
%   the single selection mode).  This will be [] when OK is 0.  OK is 1 if
%   you push the OK button, or 0 if you push the Cancel button or close the
%   figure.
%
%   Double-clicking on an item or pressing <CR> when multiple items are
%   selected has the same effect as clicking the OK button.  Pressing <CR>
%   is the same as clicking the OK button. Pressing <ESC> is the same as
%   clicking the Cancel button.
%
%   Inputs are in parameter,value pairs:
%
%   Parameter       Description
%   'ListString'    cell array of strings for the list box.
%   'SelectionMode' string; can be 'single' or 'multiple'; defaults to
%                   'multiple'.
%   'ListSize'      [width height] of listbox in pixels; defaults
%                   to [160 300].
%   'InitialValue'  vector of indices of which items of the list box
%                   are initially selected; defaults to the first item.
%   'Name'          String for the figure's title; defaults to ''.
%   'PromptString'  string matrix or cell array of strings which appears 
%                   as text above the list box; defaults to {}.
%   'OKString'      string for the OK button; defaults to 'OK'.
%   'CancelString'  string for the Cancel button; defaults to 'Cancel'.
%   'Resize'        string: can be 'on' (default) or 'off'.
%   'WindowStyle'   string: can be 'normal' or 'modal' (default).
%   'FontSize'      scalar: the Font size for texts.
%
%   A 'Select all' button is provided in the multiple selection case.
%
%   Example:
%     d = dir;
%     str = {d.name};
%     [s,v] = listdlg('PromptString','Select a file:',...
%                     'SelectionMode','single',...
%                     'ListString',str)
%
%  See also DIALOG, ERRORDLG, HELPDLG, INPUTDLG,
%    MSGBOX, QUESTDLG, WARNDLG.
%
% THIS VERSION IS MODIFIED TO USE DefaultUicontrolFontSize and be resizable

%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 1.20.4.10 $  $Date: 2009/10/24 19:19:47 $
%

%   'uh'            uicontrol button height, in pixels; default = 22.
%   'fus'           frame/uicontrol spacing, in pixels; default = 8.
%   'ffs'           frame/figure spacing, in pixels; default = 8.

% simple test:
%
% d = dir; [s,v] = listdlg('PromptString','Select a file:','ListString',{d.name});
% 

% Generate a warning in -nodisplay and -noFigureWindows mode.
% warnfiguredialog('listdlg');

error(nargchk(1,inf,nargin))

figname = mfilename;
smode = 2;   % (multiple)
promptstring = {};
liststring = [];
listsize = [160 300];
initialvalue = [];
okstring = 'OK';
cancelstring = 'Cancel';
fus = 8;
ffs = 8;
uh = 22;
resize = 'on';
windowstyle = 'modal';
fontsize = get(0,'DefaultUicontrolFontSize');
[~,tag] = fileparts(tempname); tag = [ 'ListDialogAppData_' tag ];

if mod(length(varargin),2) ~= 0
    % input args have not com in pairs, woe is me
    error('MATLAB:listdlg:InvalidArgument', 'Arguments to LISTDLG must come param/value in pairs.')
end
for i=1:2:length(varargin)
    switch lower(varargin{i})
     case 'name'
      figname = varargin{i+1};
     case 'promptstring'
      promptstring = varargin{i+1};
     case 'selectionmode'
      switch lower(varargin{i+1})
       case 'single'
        smode = 1;
       case 'multiple'
        smode = 2;
      end
     case 'listsize'
      listsize = varargin{i+1};
     case 'liststring'
      liststring = varargin{i+1};
     case 'initialvalue'
      initialvalue = varargin{i+1};
     case 'uh'
      uh = varargin{i+1};
     case 'fus'
      fus = varargin{i+1};
     case 'ffs'
      ffs = varargin{i+1};
     case 'okstring'
      okstring = varargin{i+1};
     case 'cancelstring'
      cancelstring = varargin{i+1};
     case 'resize'
      resize = varargin{i+1};
     case {'windowstyle','createmode'}
      windowstyle = varargin{i+1};
     case 'tag'
      tag = varargin{i+1};
     otherwise
      error('MATLAB:listdlg:UnknownParameter', 'Unknown parameter name passed to LISTDLG.  Name was %s', varargin{i})
    end
end

if ischar(promptstring)
    promptstring = cellstr(promptstring); 
end

if isempty(initialvalue)
    initialvalue = 1;
end

if isempty(liststring)
    error('MATLAB:listdlg:NeedParameter', 'ListString parameter is required.')
end

ex = fontsize*1.7;  % height extent per line of uicontrol text (approx)

fp = get(0,'DefaultFigurePosition');
w = 2*(fus+ffs)+listsize(1);
h = 2*ffs+6*fus+ex*length(promptstring)+listsize(2)+uh+(smode==2)*(fus+uh);
fp = [fp(1) fp(2)+fp(4)-h w h];  % keep upper left corner fixed

fig_props = { ...
    'name'                   figname ...
    'color'                  get(0,'DefaultUicontrolBackgroundColor') ...
    'resize'                 resize ...
    'numbertitle'            'off' ...
    'menubar'                'none' ...
    'windowstyle'            windowstyle ...
    'visible'                'off' ...
    'createfcn'              ''    ...
    'position'               fp   ...
    'closerequestfcn'        'delete(gcbf)', ...
    'NextPlot'               'new' ...
            };

liststring=cellstr(liststring);

fig = figure(fig_props{:});

if length(promptstring)>0
    prompt_text = uicontrol('style','text','string',promptstring,...
        'horizontalalignment','left',...
        'fontsize',fontsize, ...
        'position',[ffs+fus fp(4)-(ffs+fus+ex*length(promptstring)) ...
        listsize(1) ex*length(promptstring)]); %#ok
    set(prompt_text,'unit','normalized');
end

btn_wid = (fp(3)-2*(ffs+fus)-fus)/2;

listbox = uicontrol('style','listbox',...
                    'position',[ffs+fus ffs+uh+4*fus+(smode==2)*(fus+uh) listsize],...
                    'string',liststring,...
                    'backgroundcolor','w',...
                    'max',smode,...
                    'tag','listbox',...
                    'value',initialvalue, ...
                    'callback', {@doListboxClick});

ok_btn = uicontrol('style','pushbutton',...
                   'string',okstring,...
                   'position',[ffs+fus ffs+fus btn_wid uh],...
                   'callback',{@doOK,listbox});

cancel_btn = uicontrol('style','pushbutton',...
                       'string',cancelstring,...
                       'position',[ffs+2*fus+btn_wid ffs+fus btn_wid uh],...
                       'callback',{@doCancel,listbox});
                   
set([listbox ok_btn cancel_btn ],'Unit','normalized')

if smode == 2
    selectall_btn = uicontrol('style','pushbutton',...
                              'string','Select all',...
                              'position',[ffs+fus 4*fus+ffs+uh listsize(1) uh],...
                              'tag','selectall_btn',...
                              'callback',{@doSelectAll, listbox});
    set(selectall_btn, 'Unit','normalized');

    if length(initialvalue) == length(liststring)
        set(selectall_btn,'enable','off')
    end
    set(listbox,'callback',{@doListboxClick, selectall_btn})
end

set([fig, ok_btn, cancel_btn, listbox], 'keypressfcn', {@doKeypress, listbox});

% set(fig,'position',getnicedialoglocation(fp, get(fig,'Units')));
% Make ok_btn the default button.
% setdefaultbutton(fig, ok_btn);

% make sure we are on screen
movegui(fig)
set(fig, 'visible','on'); drawnow;

ad.value   = NaN;
ad.listbox = listbox;
ad.fig     = fig;
ad.button  = '';
ad.selection=[];
ad.selectionString='';
ad.ok_btn  = ok_btn;
ad.cancel_btn=cancel_btn;
ad.tag = tag;
set(fig, 'UserData', ad);
setappdata(0,tag, ad);

if strcmp(windowstyle,'modal')
	try
		% Give default focus to the listbox *after* the figure is made visible
		uicontrol(listbox);
		uiwait(fig);
	catch
		if ishghandle(fig)
		    delete(fig)
		end
	end

	if isappdata(0,tag)
		ad = getappdata(0,tag);
		selection = ad.selection;
		value = ad.value;
		rmappdata(0,tag)
	else
		% figure was deleted
		selection = [];
		value = 0;
	end
else % non-modal, normal
	selection = fig;
	value = ad;
end

%% figure, OK and Cancel KeyPressFcn
function doKeypress(src, evd, listbox) %#ok
switch evd.Key
 case 'escape'
  doCancel([],[],listbox);
end

%% OK callback
function doOK(ok_btn, evd, listbox) %#ok
fig = get(listbox, 'Parent');
ad = get(fig,'UserData');
if isappdata(0, ad.tag)
  ad = getappdata(0,ad.tag);
end
ad.value = 1;
ad.selection = get(listbox,'value');
string = get(listbox,'string');
ad.selectionString = string(ad.selection);
ad.button = 'OK';
set(fig,'UserData',ad);
setappdata(0,ad.tag,ad);
close(gcbf);

%% Cancel callback
function doCancel(cancel_btn, evd, listbox) %#ok
fig = get(listbox, 'Parent');
ad = get(fig,'UserData');
if isappdata(0, ad.tag)
  ad = getappdata(0,ad.tag);
end
ad.value = 0;
ad.selection = [];
string = get(listbox,'string');
ad.selectionString = string(ad.selection);
ad.button = 'Cancel';
set(fig,'UserData',ad);
setappdata(0,ad.tag,ad)
close(gcbf);

%% SelectAll callback
function doSelectAll(selectall_btn, evd, listbox) %#ok
set(selectall_btn,'enable','off')
set(listbox,'value',1:length(get(listbox,'string')));

%% Listbox callback
function doListboxClick(listbox, evd, selectall_btn) %#ok
% if this is a doubleclick, doOK
if strcmp(get(gcbf,'SelectionType'),'open')
    doOK([],[],listbox);
elseif nargin == 3
    if length(get(listbox,'string'))==length(get(listbox,'value'))
        set(selectall_btn,'enable','off')
    else
        set(selectall_btn,'enable','on')
    end
end
