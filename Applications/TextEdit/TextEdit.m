function hF = TextEdit(filename, options)
% Simple Text editor
%
%  TextEdit
%    Open an empty editor window
%  TextEdit(filename)
%    Open the file content inside a new editor window
%  TextEdit(string)
%    Display the string inside a new Editor window
%  TextEdit(..., name)
%    Give the TextEdit a name
%  TextEdit(filename|string, options)
%    with options as a structure:
%      name:    the name of the TextEdit
%      command: when 1, displays a command line, in which case the main pane is not
%               editable. Commands entered in this area are echoed in the main 
%               pane and executed.
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
if ~isfield(options, 'command') || isempty(options.command)
  options.command = 0;
end

% build the main window
hF=figure('MenuBar','none',...
    'Name',options.name,'Resize','on',...
    'Position',[0 0 600 400],'Color','w','CloseRequestFcn', @textedit_quit);
centerfig();

% Menu File
hMA=uimenu(hF,'Label','File');
uimenu(hMA,'Label','New','Callback','TextEdit', 'Accelerator','n');
uimenu(hMA,'Label','Open...','Callback',@textedit_open, 'Accelerator','o');
uimenu(hMA,'Label','Save...','Callback',@textedit_save, 'Accelerator','s');
uimenu(hMA,'Label','Evaluate','Callback',@textedit_eval, 'Accelerator','e');
uimenu(hMA,'Label','Quit','Callback',@textedit_quit,'Separator','on','Accelerator','q');

% Menu Edit
hME=uimenu(hF,'Label','Edit');
uimenu(hME,'Label','Cut','Callback',@textedit_cut);
uimenu(hME,'Label','Copy','Callback',@textedit_copy);
uimenu(hME,'Label','Paste','Callback',@textedit_paste);

% Menu Help
hMA=uimenu(hF,'Label','Help');
uimenu(hMA,'Label','About...','Callback',@textedit_about, 'Accelerator','h');
    
use_fallback = 1;
% we use the nice Java SyntaxTextPane
try
  if usejava('jvm') && usejava('swing') && feature('ShowFigureWindows')
    jCodePane = com.mathworks.widgets.SyntaxTextPane;
    codeType = jCodePane.M_MIME_TYPE;  % jCodePane.contentType='text/m-MATLAB'
    jCodePane.setContentType(codeType)
    jCodePane.setText('')
    jScrollPane = com.mathworks.mwswing.MJScrollPane(jCodePane);
    [jhPanel,hContainer] = javacomponent(jScrollPane,[10,10,300,100],hF);
    hTxt = jCodePane;
    % set the editor panel to whole figure, normalised.
    set(hContainer,'Units',   'normalized')
 
    if options.command
      set(hContainer,'Position',[0 0.1 1 0.9])
      % we create an other ScrollPane with highlight syntax and a 'Run' button
      jCodePane2 = com.mathworks.widgets.SyntaxTextPane;
      codeType = jCodePane2.M_MIME_TYPE;  % jCodePane.contentType='text/m-MATLAB'
      jCodePane2.setContentType(codeType)
      str = ['% enter Matlab commands here\n' ];
      str = sprintf(strrep(str,'%','%%'));
      jCodePane2.setText(str)
      jScrollPane2 = com.mathworks.mwswing.MJScrollPane(jCodePane2);
      [jhPanel2,hContainer2] = javacomponent(jScrollPane2,[10,10,300,100],hF);
      set(hContainer2,'Units',   'normalized')
      set(hContainer2,'Position',[0 0 1 0.1]);
      options.command_pane_java = jCodePane2;
      options.command_pane_uicontrol = hContainer2;
    else
      set(hContainer,'Position',[0 0 1 1])
      options.command_pane = [];
      options.command_pane_uicontrol = [];
    end
    options.display_pane = jCodePane;
    options.display_pane_uicontrol = hContainer;
    
    % more menu stuff
    uimenu(hME,'Label','Select All', 'Callback','tmp_hTxt=get(gcbf,''UserData''); tmp_hTxt.display_pane.selectAll; clear tmp_hTxt;')
    use_fallback = 0;
  end
catch
  use_fallback = 1;
end
if use_fallback
  % fallback Editor 
  hTxt=uicontrol('style','edit','String','',...
    'Units','Normalized','Position',[0 0 1 1],...
    'BackgroundColor','w','Horizontal','left',...
    'Tag','TextEdit_Text', ...
    'Max',1000,'FontSize',10,'FontName','Arial');
    
  % Align menu
  h_AT=uimenu(hME,'Label','Align');
  uimenu(h_AT,'Label','Left','Callback',@textedit_align);
  uimenu(h_AT,'Label','Right','Callback',@textedit_align);
  uimenu(h_AT,'Label','Center','Callback',@textedit_align);

  % Theme selection
  hMT=uimenu(hF,'Label','Layout');
  uimenu(hMT,'Label','Blue','Callback',@textedit_theme);
  uimenu(hMT,'Label','Classic','Callback',@textedit_theme);
  uimenu(hMT,'Label','Cool','Callback',@textedit_theme);
  uimenu(hMT,'Label','Dark','Callback',@textedit_theme);
  uimenu(hMT,'Label','Silver','Callback',@textedit_theme);
  
  % more in Edit menu
  uimenu(hME,'Label','Change Font Color...','Callback',@textedit_fontcolor, 'Separator','on');
  uimenu(hME,'Label','Change Font...','Callback',@textedit_font);
  options.display_pane = hTxt;
  options.command = 0;
end

% protect figure from over-plotting
set(hF,'HandleVisibility','callback','UserData',options);

if nargin 
  if iscellstr(filename)
    filename = char(filename);
  elseif isa(filename, 'function_handle')
    filename = func2str(filename);
  end
end

if nargin && ischar(filename)
  if size(filename,1) == 1 && ~isempty(dir(filename))
    textedit_load(filename);
  else
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

% Save content into a text file
    function textedit_save(~,~)
        txt = textedit_getText(options.display_pane);
        [filename, pathname] = uiputfile( ...
          {'*.txt','Text files (*.txt)'; ...
           '*.m','Matlab script (*.m)'; ...
           '*.*','All files'}, 'TextEdit: Save as:');
        if isequal(filename,0) || isequal(pathname,0)
            return;
        else
            fid=fopen(fullfile(pathname,filename),'wt');
            if fid == -1
              error([ mfilename ': Could not open file ' fullfile(pathname,filename) ]);
            end
            for i=1:size(txt,1)
                fprintf(fid,'%s\n',txt(i,:));
            end
            fclose(fid);
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
        disp([ '% ' mfilename ': Evaluating code: ' datestr(now) ])
        disp(txt)
        disp([ '% ' mfilename ': end of code to evaluate.' ])
        disp(' ')
        try
          evalin('base', txt');
        catch ME
          disp(ME.message)
        end
    end
    
    function textedit_quit(~,~)
      options = get(hF,'UserData');
      try
        delete(options.display_pane);
      end
      delete(hF);
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

% Select font color
    function textedit_fontcolor(~,~)
      if ~isa(options.display_pane, 'com.mathworks.widgets.SyntaxTextPane')
        clr=uisetcolor([1 1 1]);
        set(options.display_pane,'ForegroundColor',clr);
      end
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

% Align text
    function textedit_align(src,~)
      if ~isa(options.display_pane, 'com.mathworks.widgets.SyntaxTextPane')
        alin=get(src,'Label');
        set(options.display_pane,'Horizontal',alin);
      end
    end

% Modify the Theme/color set
    function textedit_theme(src,~)
        tipo=get(src,'Label');
        MTEM={'Blue',[0.2 0.2 0.95],[1 1 1];
            'Classic',[1 1 1],[0 0 1];
            'Cool',[0.9 0.5 0.9],[1 1 0];
            'Dark',[0 0 0],[1 1 1];
            'Silver',[0.6 0.6 0.6],[1 1 1]};
        k=strfind(MTEM(:,1),tipo);
        k=cellfun(@isempty,k);
        if ~isa(options.display_pane, 'com.mathworks.widgets.SyntaxTextPane')
          set(options.display_pane,'BackgroundColor',MTEM{find(~k),2});
          set(options.display_pane,'ForegroundColor',MTEM{find(~k),3});
        end
    end

% About
    function textedit_about(~,~)

        devel='By: Jorge De Los Santos; Improved by E. Farhi, ILL, 2017.';
        e_mail='E-mail: delossantosmfq@gmail.com; farhi@ill.fr';
        blog='Blog: http://matlab-typ.blogspot.mx ; http://ifit.mccode.org';
        nvrs='TextEdit 1.0';
        
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
