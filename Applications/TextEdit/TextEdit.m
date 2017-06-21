function hF = TextEdit(filename, name)
% Simple Text editor
%
%  TextEdit
%    Open an empty editor window
%  TextEdit(filename)
%    Open the file content inside a new editor window
% TextEdit(string)
%    Display the string inside a new Editor window
%
% $ Por: Jorge De Los Santos $
% $ E-mail: delossantosmfq@gmail.com $
% $ Blog: http://matlab-typ.blogspot.mx $
% $ Rev. 0.0.1 $ 02/08/2014 $
%
% http://fr.mathworks.com/matlabcentral/fileexchange/47614-textedit
% Modified by E. Farhi, 2016 into english

if nargin < 2, name=mfilename; end
hF=figure('MenuBar','none',...
    'Name',name,'Resize','on',...
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
    
% we use the nice Java SyntaxTextPane
try
  jCodePane = com.mathworks.widgets.SyntaxTextPane;
  codeType = jCodePane.M_MIME_TYPE;  % jCodePane.contentType='text/m-MATLAB'
  jCodePane.setContentType(codeType)
  str = ['% enter Matlab code here\n' ];
  str = sprintf(strrep(str,'%','%%'));
  jCodePane.setText(str)
  jScrollPane = com.mathworks.mwswing.MJScrollPane(jCodePane);
  [jhPanel,hContainer] = javacomponent(jScrollPane,[10,10,300,100],hF);
  hTxt = jCodePane;
  % set the editor panel to whole figure, normalised.
  set(hContainer,'Units',   'normalized')
  set(hContainer,'Position',[0 0 1 1])
  
  % more menu stuff
  uimenu(hME,'Label','Select All', 'Callback','tmp_hTxt=get(gcbf,''UserData''); tmp_hTxt.selectAll; clear tmp_hTxt;')
catch
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
end

% protect figure from over-plotting
set(hF,'HandleVisibility','callback','UserData',hTxt);

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
    textedit_setText(hTxt, filename)
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
        textedit_setText(hTxt, txt);
        if isa(hTxt, 'com.mathworks.widgets.SyntaxTextPane')
          hTxt.setToolTipText(sprintf('File: %s\nSize: %i', ...
            filename, numel(txt)));
        elseif ishandle(hTxt)
          set(hTxt, 'ToolTip', sprintf('File: %s\nSize: %i', ...
            filename, numel(txt)));
        end
        set(hF,'Name', [ 'TextEdit: ' filename ]);
    end

% Save content into a text file
    function textedit_save(~,~)
        txt = textedit_getText(hTxt);
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
        if isa(hTxt,'com.mathworks.widgets.SyntaxTextPane')
          txt = char(hTxt.getSelectedText);
        else
          txt = [];
        end
        if isempty(txt)
          txt = char(textedit_getText(hTxt));
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
      hTxt = get(hF,'UserData');
      try
        delete(hTxt);
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
        hTxt = get(gcbf,'UserData');
        clipboard('copy', textedit_getText(hTxt));
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
        hTxt = get(gcbf,'UserData');
        str1 = textedit_getText(hTxt);
        str2 = clipboard('paste');
        textedit_setText(hTxt, strvcat(str1, str2));
      end
    end

%  Copy and remove the whole content
    function textedit_cut(~,~)
        str = textedit_getText(hTxt);
        clipboard('copy', str);
        textedit_setText(hTxt, '');
    end

% Select font color
    function textedit_fontcolor(~,~)
      if ~isa(hTxt, 'com.mathworks.widgets.SyntaxTextPane')
        clr=uisetcolor([1 1 1]);
        set(hTxt,'ForegroundColor',clr);
      end
    end

% Select Font
    function textedit_font(~,~)
        if ~isa(hTxt, 'com.mathworks.widgets.SyntaxTextPane')
          fmod=uisetfont();
          if ~isequal(fmod,0)
            set(hTxt,'FontName',fmod.FontName,...
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
      if ~isa(hTxt, 'com.mathworks.widgets.SyntaxTextPane')
        alin=get(src,'Label');
        set(hTxt,'Horizontal',alin);
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
        if ~isa(hTxt, 'com.mathworks.widgets.SyntaxTextPane')
          set(hTxt,'BackgroundColor',MTEM{find(~k),2});
          set(hTxt,'ForegroundColor',MTEM{find(~k),3});
        end
    end

% About
    function textedit_about(~,~)
        figure('MenuBar','none','NumberTitle','off',...
            'Name','TextEdit: About','Resize','off',...
            'Position',[0 0 200 100],'color','w');
        centerfig();
        devel='By: Jorge De Los Santos';
        e_mail='E-mail: delossantosmfq@gmail.com';
        blog='Blog: http://matlab-typ.blogspot.mx';
        nvrs='TextEdit 0.0.2';
        uicontrol('style','text','String',devel,...
            'Units','Normalized','Position',[0.1 0.80 0.8 0.15],...
            'FontName','Arial Narrow','FontSize',10,...
            'ForegroundColor',ones(1,3)*0.2);
        uicontrol('style','text','String',{e_mail,blog},...
            'Units','Normalized','Position',[0.1 0.45 0.8 0.3],...
            'FontName','Arial Narrow','FontSize',9,...
            'ForegroundColor',ones(1,3)/2);
        uicontrol('style','text','String',nvrs,...
            'Units','Normalized','Position',[0.1 0.15 0.8 0.1],...
            'FontName','Courier','FontSize',10,'FontWeight','b',...
            'ForegroundColor',[0 0 0.5]);
        set(findobj('style','text'),'BackgroundColor','w');
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
