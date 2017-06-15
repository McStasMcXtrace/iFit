function hF = TextEdit(filename)
% Simpl Text editor
% 
%
% $ Por: Jorge De Los Santos $
% $ E-mail: delossantosmfq@gmail.com $
% $ Blog: http://matlab-typ.blogspot.mx $
% $ Rev. 0.0.1 $ 02/08/2014 $
%
% http://fr.mathworks.com/matlabcentral/fileexchange/47614-textedit
% Modified by E. Farhi, 2016 into english

hF=figure('MenuBar','none',...
    'Name','TextEdit','Resize','on',...
    'Position',[0 0 600 400],'Color','w');
centerfig();

% Menu File
hMA=uimenu(hF,'Label','File');
uimenu(hMA,'Label','New','Callback','TextEdit', 'Accelerator','n');
uimenu(hMA,'Label','Open...','Callback',@textedit_open, 'Accelerator','o');
uimenu(hMA,'Label','Save...','Callback',@textedit_save, 'Accelerator','s');
uimenu(hMA,'Label','Evaluate','Callback',@textedit_eval, 'Accelerator','e');
uimenu(hMA,'Label','Quit','Callback','delete(gcbf)','Separator','on','Accelerator','q');

% Menu Edit
hME=uimenu(hF,'Label','Edit');
uimenu(hME,'Label','Cut','Callback',@textedit_cut);
uimenu(hME,'Label','Copy','Callback',@textedit_copy);
uimenu(hME,'Label','Paste','Callback',@textedit_paste);
uimenu(hME,'Label','Change Font Color...','Callback',@textedit_fontcolor, 'Separator','on');
uimenu(hME,'Label','Change Font...','Callback',@textedit_font);
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

% Menu Help
hMA=uimenu(hF,'Label','Help');
uimenu(hMA,'Label','About...','Callback',@textedit_about, 'Accelerator','h');

% Editor 
hTxt=uicontrol('style','edit','String','',...
    'Units','Normalized','Position',[0 0 1 1],...
    'BackgroundColor','w','Horizontal','left',...
    'Tag','TextEdit_Text', ...
    'Max',1000,'FontSize',10,'FontName','Arial');

% protect figure from over-plotting
set(hF,'HandleVisibility','callback','UserData',hTxt);

if nargin && ischar(filename)
  textedit_load(filename);
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
        set(hTxt, 'String',txt);
        set(hTxt, 'ToolTip', sprintf('File: %s\nSize: %i', ...
          filename, numel(txt)));
        set(hF,'Name', [ 'TextEdit: ' filename ]);
    end

% Save content into a text file
    function textedit_save(~,~)
        txt=get(hTxt,'String');
        [filename, pathname] = uiputfile({'*.txt','Text files (*.txt)'; ...
        '*.*','All files'}, 'TextEdit: Save text as:');
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
        txt=get(hTxt,'String');
        try
          evalin('base', txt');
        catch ME
          disp(ME.message)
        end
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
        clipboard('copy', get(hTxt,'String'));
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
        str1 = get(hTxt,'String');
        str2 = clipboard('paste');
        set(hTxt,'String',strvcat(str1, str2));
      end
    end

%  Copy and remove the whole content
    function textedit_cut(~,~)
        str = get(hTxt,'String');
        clipboard('copy', str);
        set(hTxt,'String','');
    end

% Select font color
    function textedit_fontcolor(~,~)
        clr=uisetcolor([1 1 1]);
        set(hTxt,'ForegroundColor',clr);
    end

% Select Font
    function textedit_font(~,~)
        fmod=uisetfont();
        if ~isequal(fmod,0)
            set(hTxt,'FontName',fmod.FontName,...
                'FontSize',fmod.FontSize,...
                'FontWeight',fmod.FontWeight,...
                'FontAngle',fmod.FontAngle);
        else
            return;
        end
    end

% Align text
    function textedit_align(src,~)
        alin=get(src,'Label');
        set(hTxt,'Horizontal',alin);
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
        set(hTxt,'BackgroundColor',MTEM{find(~k),2});
        set(hTxt,'ForegroundColor',MTEM{find(~k),3});
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
