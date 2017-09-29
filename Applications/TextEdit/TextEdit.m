function hF = TextEdit(filename, options)
% Simple Text editor
%
%  TextEdit
%    Open an empty editor window
%  TextEdit(filename)
%    Open the file content inside a new editor window
%  TextEdit(string or cellstr)
%    Display the string inside a new Editor window
%  TextEdit(..., name)
%    Give the TextEdit a name
%
% $ Initially from: Jorge De Los Santos $
% $ E-mail: delossantosmfq@gmail.com $
% $ Blog: http://matlab-typ.blogspot.mx $
% $ Rev. 0.0.1 $ 02/08/2014 $
%
% Improved with Java component thanks to Yair Altman comments.
% by E. Farhi <farhi@ill.fr> Institut Laue-Langevin BSD, 2017
% Rev 1.0
%
% Orig: http://fr.mathworks.com/matlabcentral/fileexchange/47614-textedit
% Modified by E. Farhi, 2016 into english

% build the 'options'
if nargin < 2, options.name=[]; end
if ischar(options) name = options; options=''; options.name = name; end

if ~isstruct(options) || ~isfield(options, 'name') || isempty(options.name), 
  options.name    = mfilename;
end

% build the main window
hF=figure('MenuBar','none',...
    'Tag','TextEdit', ...
    'Name',options.name,'Resize','on',...
    'Position',[0 0 600 400],'Color','w','CloseRequestFcn', 'delete(gcbf)');
centerfig();

% Menu File
hMA=uimenu(hF,'Label','File');
uimenu(hMA,'Label','New','Callback','TextEdit', 'Accelerator','n');
uimenu(hMA,'Label','Open...','Callback',@textedit_open, 'Accelerator','o');
uimenu(hMA,'Label','Save','Callback',@textedit_save);
uimenu(hMA,'Label','Save as...','Callback',@textedit_saveas, 'Accelerator','s');
uimenu(hMA,'Label','Evaluate selection','Callback',@textedit_eval, 'Accelerator','e');

% Menu Edit
hME=uimenu(hF,'Label','Edit');
uimenu(hME,'Label','Cut','Callback',@textedit_cut);
uimenu(hME,'Label','Copy','Callback',@textedit_copy);
uimenu(hME,'Label','Paste','Callback',@textedit_paste);

% Menu Help
hMH=uimenu(hF,'Label','Help');
if isdeployed
  uimenu(hMH,'Label','iFit Terminal','Callback','doc(iData,''iFit'')');
end
uimenu(hMH,'Label','Main iFit help','Callback','doc(iData,''index'')');
uimenu(hMH,'Label','Data set object (iData)','Callback','doc(iData,''iData'')');
uimenu(hMH,'Label','Data set methods (iData)','Callback','doc(iData,''Methods'')');
uimenu(hMH,'Label','Data set Math','Callback','doc(iData,''Math'')');
uimenu(hMH,'Label','Model object (iFunc)','Callback','doc(iData,''iFunc'')');
uimenu(hMH,'Label','Predefined Models (iFunc)','Callback','doc(iData,''Models'')');
uimenu(hMH,'Label','Importing data','Callback','doc(iData,''Load'')');
uimenu(hMH,'Label','Exporting data','Callback','doc(iData,''Save'')');
uimenu(hMH,'Label','Plotting data','Callback','doc(iData,''Plot'')');
uimenu(hMH,'Label','Fitting data/model','Callback','doc(iData,''Fit'')');
uimenu(hMH,'Label','About...','Callback',@textedit_about, 'Accelerator','h','separator','on');

% add a toolbar
iconsroot = fullfile(matlabroot,'toolbox','matlab','icons');
tools = { ...
  'New',  'file_new.png', 'TextEdit'; ...
  'Open', 'file_open.png', @textedit_open; ...
  'Save', 'file_save.png', @textedit_save; ...
  'Run',  'greenarrowicon.gif', @textedit_eval; ...
  'Help', 'helpicon.gif',  'doc(iData,''index'')' };
ht=uitoolbar(hF);
for index=1:size(tools,1)
  [X map] = imread(fullfile(iconsroot, [ tools{index,2} ]));
  if ~isempty(map)
    icon = ind2rgb(X,map);
  elseif ndims(X) == 3 && size(X,3) == 3
    maxX = double(max(X(:)));
    icon = double(X)/maxX;
  end
  icon(icon == 0) = nan;
  uipushtool(ht,'CData', icon, ...
                 'TooltipString',tools{index,1},...
                 'ClickedCallback',tools{index,3});
end
    
use_fallback = 1;
% we use the nice Java SyntaxTextPane
try
  if usejava('jvm') && usejava('swing') && feature('ShowFigureWindows')
    jCodePane = com.mathworks.widgets.SyntaxTextPane;
    codeType = jCodePane.M_MIME_TYPE;  % jCodePane.contentType='text/m-MATLAB'
    jCodePane.setContentType(codeType)
    jCodePane.setText('')
    jCodePane.setDragEnabled(true)
    jCodePane.setTextDragAndDropEnabled(true)
    jScrollPane = com.mathworks.mwswing.MJScrollPane(jCodePane);
    [jhPanel,hContainer] = javacomponent(jScrollPane,[10,10,300,100],hF);
    hTxt = jCodePane;
    % set the editor panel to whole figure, normalised.
    set(hContainer,'Units',   'normalized')
    set(hContainer,'Position',[0 0 1 1])

    options.display_pane = jCodePane;
    options.display_pane_uicontrol = hContainer;
    
    % more menu stuff
    uimenu(hME,'Label','Select All', 'Callback', ...
      'tmp_hTxt=get(gcbf,''UserData''); tmp_hTxt.display_pane.selectAll; clear tmp_hTxt;')
    % print is broken: blocks the jCodePane
    % uimenu(hMA,'Label','Print','Callback', ...
    %  'tmp_hTxt=get(gcbf,''UserData''); tmp_hTxt.display_pane.print; clear tmp_hTxt;', ...
    %  'Separator','on','Accelerator','p');
    use_fallback = 0;
  end
catch ME
  disp(getReport(ME))
  use_fallback = 1;
end
if use_fallback
  % fallback Editor 
  hTxt=uicontrol('style','edit','String','',...
    'Units','Normalized','Position',[0 0 1 1],...
    'BackgroundColor','w','Horizontal','left',...
    'Tag','TextEdit_Text', ...
    'Max',1000,'FontSize',10,'FontName','Arial');
  
  % more in Edit menu
  uimenu(hME,'Label','Change Font...','Callback',@textedit_font);
  options.display_pane = hTxt;
end

% add Quit
uimenu(hMA,'Label','Quit','Callback','delete(gcbf)','Separator','on','Accelerator','q');

% protect figure from over-plotting
set(hF,'HandleVisibility','callback','UserData',options);

if nargin 
  if iscellstr(filename)
    filename = char(filename);
  elseif isa(filename, 'function_handle')
    filename = func2str(filename);
  end
end

options.filename = '';
if nargin && ischar(filename)
  try
      if size(filename,1) == 1 && ~isempty(dir(filename))
        textedit_load(filename);
        options.filename = filename;
      else
        textedit_setText(options.display_pane, filename)
      end
  catch
      textedit_setText(options.display_pane, filename)
  end
end

% - - - - - - - - - - - - - Functions - - - - - - - - - - - - - - - - - - 

% Open a file
    function textedit_open(~,~)
        [filename,pathname] = uigetfile({'*.txt','Text files (*.txt)'; ...
        '*.*','All files'}, 'TextEdit: Select a text file:');
        if isequal(filename,0) || isequal(pathname,0)
            return;
        else
            textedit_load(fullfile(pathname,filename));
        end
    end
    
    function textedit_load(filename)
        txt=fileread(filename);
        textedit_setText(options.display_pane, txt);
        if isa(options.display_pane, 'com.mathworks.widgets.SyntaxTextPane')
          options.display_pane.setToolTipText(sprintf('File: %s\nSize: %i', ...
            filename, numel(txt)));
        elseif ishandle(options.display_pane)
          set(options.display_pane, 'ToolTip', sprintf('File: %s\nSize: %i', ...
            filename, numel(txt)));
        end
        set(hF,'Name', [ 'TextEdit: ' filename ]);
    end

% Save content into an existing text file
    function textedit_save(~,~)
      if isempty(options.filename) || ~ischar(options.filename)
        textedit_saveas;
        return;
      else
        txt = textedit_getText(options.display_pane);
        fid=fopen(options.filename,'wt');
        if fid == -1
          error([ mfilename ': Could not open file ' options.filename ' for saving.' ]);
        end
        for i=1:size(txt,1)
            fprintf(fid,'%s\n',char(txt(i,:)));
        end
        fclose(fid);
        set(hF,'Name', [ 'TextEdit: ' options.filename ]);
      end
    end
    
% Save content into a (new) text file
    function textedit_saveas(~,~)
        txt = textedit_getText(options.display_pane);
        if ~isempty(options.filename) && ischar(options.filename)
          [filename, pathname] = uiputfile( ...
            {'*.txt','Text files (*.txt)'; ...
             '*.m','Matlab script (*.m)'; ...
             '*.*','All files'}, 'TextEdit: Save as:', options.filename);
        else
          [filename, pathname] = uiputfile( ...
            {'*.txt','Text files (*.txt)'; ...
             '*.m','Matlab script (*.m)'; ...
             '*.*','All files'}, 'TextEdit: Save as:');
        end
        if isequal(filename,0) || isequal(pathname,0)
            return;
        else
            fid=fopen(fullfile(pathname,filename),'wt');
            if fid == -1
              error([ mfilename ': Could not open file ' fullfile(pathname,filename) ' for saving.' ]);
            end
            for i=1:size(txt,1)
                fprintf(fid,'%s\n',char(txt(i,:)));
            end
            fclose(fid);
            set(hF,'Name', [ 'TextEdit: ' fullfile(pathname,filename) ]);
            options.filename = fullfile(pathname,filename);
        end
    end
    
% Evaluate whole text
    function textedit_eval(~,~)
        % get Selected text
        if isa(options.display_pane,'com.mathworks.widgets.SyntaxTextPane')
          txt = char(options.display_pane.getSelectedText);
        else
          txt = [];
        end
        if isempty(txt)
          txt = char(textedit_getText(options.display_pane));
        end
        if isempty(txt), return; end
        disp([ '% ' mfilename ': Evaluating code: ' datestr(now) ' from Figure ' num2str(get(hF,'Name')) ])
        disp(txt)
        disp([ '% ' mfilename ': end of code to evaluate.' ])
        disp(' ')
        set(hF, 'Pointer', 'watch')
        try
          evalin('base', txt');
        catch ME
          disp(ME.message)
        end
        set(hF, 'Pointer', 'arrow')
    end

% Copy into the clipboard
    function textedit_copy(~,~)
      try
        % Use java.awt.Robot to simulate the Ctrl-C/V keys
        import java.awt.Robot;
        import java.awt.event.KeyEvent;
        rb=Robot();
        rb.keyPress(KeyEvent.VK_CONTROL);
        rb.keyPress(KeyEvent.VK_C);
        rb.keyRelease(KeyEvent.VK_CONTROL);
        rb.keyRelease(KeyEvent.VK_C);
        clear('rb');
      catch
        options = get(gcbf,'UserData');
        clipboard('copy', textedit_getText(options.display_pane));
      end
    end

% Paste from the clipboard
    function textedit_paste(~,~)
      try
        % Use java.awt.Robot to simulate the Ctrl-C/V keys
        import java.awt.Robot;
        import java.awt.event.KeyEvent;
        rb=Robot();
        rb.keyPress(KeyEvent.VK_CONTROL);
        rb.keyPress(KeyEvent.VK_V);
        rb.keyRelease(KeyEvent.VK_CONTROL);
        rb.keyRelease(KeyEvent.VK_V);
        clear('rb');
      catch
        options = get(gcbf,'UserData');
        str1 = textedit_getText(options.display_pane);
        str2 = clipboard('paste');
        textedit_setText(options.display_pane, strvcat(str1, str2));
      end
    end

%  Copy and remove the whole content
    function textedit_cut(~,~)
        options = get(gcbf,'UserData');
        str = textedit_getText(options.display_pane);
        clipboard('copy', str);
        textedit_setText(options.display_pane, '');
    end

% Select Font
    function textedit_font(~,~)
        if ~isa(options.display_pane, 'com.mathworks.widgets.SyntaxTextPane')
          fmod=uisetfont();
          if ~isequal(fmod,0)
            set(options.display_pane,'FontName',fmod.FontName,...
                'FontSize',fmod.FontSize,...
                'FontWeight',fmod.FontWeight,...
                'FontAngle',fmod.FontAngle);
          end
        else
            return;
        end
    end


% About
    function textedit_about(~,~)

        devel='Initial design by By: Jorge De Los Santos; Improved by E. Farhi, ILL, 2017.';
        e_mail='E-mail: delossantosmfq@gmail.com; farhi@ill.fr';
        blog='Blog: http://matlab-typ.blogspot.mx ; http://ifit.mccode.org';
        nvrs='TextEdit 2.0';
        
        str = { devel, e_mail, blog, nvrs };
        helpdlg(str, nvrs);
    end
end

% set/get text from both Java or Matlab-standard uicontrol
function txt=textedit_getText(hTxt)
  if isa(hTxt,'com.mathworks.widgets.SyntaxTextPane')
    txt=hTxt.getText;
  elseif ishandle(hTxt)
    txt=get(hTxt,'String');
  end
end

function textedit_setText(hTxt, str)
  if isa(hTxt,'com.mathworks.widgets.SyntaxTextPane')
    
    if ~iscell(str)
      if isnumeric(str), str=num2str(str); end
      str = char(str);
      if size(str,1) > 1, 
        str = cellstr(str);
      end
    end
    if iscellstr(str), str=sprintf('%s\n', str{:}); end
    hTxt.setText(str);
  elseif ishandle(hTxt)
    set(hTxt, 'String',str);
  end
end
