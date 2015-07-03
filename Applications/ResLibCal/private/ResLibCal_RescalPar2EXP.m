function [EXP, titl] = ResLibCal_RescalPar2EXP(str, EXP)
% [EXP = ResLibCal_RescalPar2EXP(str): Convert a structure into a ResLib EXP
% 
% searches for ResCal parameter fields from structrue 'str', and fill in a ResLib EXP
% 
% Returns:
%  EXP: ResLib structure

% Calls: str2struct, ResLibCal_EXP2RescalPar

persistent labels

if nargin < 1, return; end
if nargin < 2, EXP = []; end
if ischar(str)
  % remove any 'ResCal:' keyword
  index = strfind(lower(str), 'rescal:');
  str(index:(index+6))=[];
  index = strfind(lower(str), 'rescal');
  str(index:(index+5))=[];
  content = str2struct(str);  % do we directly import an EXP from char ?
else 
  content = str;
end
titl = [];

if isfield(content, 'EXP')    % imported a full ResLibCal structure
  EXP = content.EXP;
  str = '';
end
if isfield(content, 'sample') % imported an EXP structure
  EXP = content;
  str = '';
end

% signification of ResCal fields
  if isempty(labels) % only the first time, then stored as persistent
    [p, labels] = ResLibCal_EXP2RescalPar([]); % get ResCal5 field names
    labels = strtok(labels);
    % add some aliases
    labels{end+1} = 'H';
    labels{end+1} = 'K';
    labels{end+1} = 'L';
    labels{end+1} = 'W';
  end
  
  % is this a ResCal structure ?
  if any(isfield(content, labels)), str = content; end
  if isempty(str), return; end

  % convert string to structure
  if ischar(str)
    % we search for <tokens>=<value> in the string
    [start,ending,match]=regexp(str, ...
      strcat('\<',labels,'\>(\s)*='),'start','end','match','once');
    p = [];
    for index=find(~cellfun(@isempty, start))
      len   = min(ending{index}+20, length(str));
      value = strtrim(str((ending{index}+1):len)); % remove blanks, pass '='
      value = str2double(strtok(value,[' ;:,' sprintf('\n\t\f') ])); % extract next word
      if isfinite(value)
        p.(labels{index}) = value;
      end
    end
    if ~isempty(p)
      str = p; % now a structure
    end
  end

  if ischar(str) && ~isempty(strfind(str, 'Title (max.60 characters)'))
    % legacy ResTrax configuration file. Very partil set of parameters.
    lines = textscan(str, '%s', 'delimiter', '\n'); % split all lines
    lines = lines{1};
    ResTrax.Source = str2double(regexp(lines{4}, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match'));
    ResTrax.Monok  = str2double(regexp(lines{8}, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match'));
    ResTrax.Ana    = str2double(regexp(lines{10}, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match'));
    ResTrax.Det    = str2double(regexp(lines{12}, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match'));
    ResTrax.Dist   = str2double(regexp(lines{14}, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match'));
    ResTrax.C1     = str2double(regexp(lines{16}, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match'));
    ResTrax.C2     = str2double(regexp(lines{18}, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match'));
    ResTrax.C3     = str2double(regexp(lines{20}, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match'));
    ResTrax.C4     = str2double(regexp(lines{22}, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match'));
    ResTrax.Collimators = str2double(regexp(lines{23}, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', 'match'));
    % now import data specifically as a ResCal structure
    if ResTrax.Source(1) == 0
      ResTrax.Source(3:4) = ResTrax.Source(2)*1.41;
    end
    p.WB=ResTrax.Source(3); p.HB=ResTrax.Source(4);
    p.TM=ResTrax.Monok(4);  p.WM=ResTrax.Monok(5); p.HM=ResTrax.Monok(6);
    p.TA=ResTrax.Ana(4);    p.WA=ResTrax.Ana(5);   p.HA=ResTrax.Ana(6);
    if ResTrax.Det(1) == 0
      ResTrax.Det(3:4) = ResTrax.Det(2)*1.41;
    end
    p.WD=ResTrax.Det(3); p.HD=ResTrax.Det(4);
    p.L1=ResTrax.Dist(1); p.L2=ResTrax.Dist(2); p.L3=ResTrax.Dist(3); p.L4=ResTrax.Dist(4);
    titl = [ 'ResTrax legacy configuration: ' strtrim(lines{2}) ];
  end

  % handle input as a numerical vector: ResCal file
  if isnumeric(str) && isvector(str)
    if numel(str) >= 42 % legacy ResCal5 .par file
      if  numel(str) > 42
          disp([ mfilename ': initial parameter file contains ' ...
              num2str(numel(str)) '. Only using up to 42.']);
          str = str(1:42);
      end
      labs=labels(1:42);
      str = mat2cell(str(:),ones(1,length(str)));
      str = cell2struct(str(:),labs(:),1);  % make a structure
      titl = 'ResCal Cooper-Nathans parameters';
    elseif numel(str) == 27
      labs=labels(42+(1:27));
      str = mat2cell(str(:),ones(1,length(str)));
      str = cell2struct(str(:),labs(:),1);  % make a structure
      titl = 'ResCal Popovici parameters (Rescal5)';
    end
  end
  
  % field names to search for are in 'labels'

  if isstruct(str)
    % field names provided in input structure  
    fields = fieldnames(str);
    p=[];
    for index=1:length(fields)
      found=find(~cellfun(@isempty, strfind(labels, fields{index})));
      if length(found) > 1
        % check if one of the matches is exact
        if length(find(strcmp(fields{index}, labels(found))))==1
          found=found(strcmp(fields{index}, labels(found)));
        else
          disp([ mfilename ': Token ' fields{index} ' matches more than one ResCal parameter.'])
          disp(labels(found))
          disp([ 'Using first match ' labels{found(1)} ])
          found = found(1);
        end
      elseif isempty(found), continue; end
      p.(labels{found}) = str.(fields{index});
    end
    if isempty(titl), titl = 'Restrax/Rescal parameters'; end
  end
  % now convert ResCal clean 'p' structure to ResLib EXP

  % ResCal parameters (pres: 42)
  if isfield(p,'DM'),   EXP.mono.d     = p.DM; end
  if isfield(p,'DA'),   EXP.ana.d      = p.DA; end
  if isfield(p,'ETAM'), EXP.mono.mosaic= p.ETAM; end
  if isfield(p,'ETAA'), EXP.ana.mosaic = p.ETAA; end
  if isfield(p,'ETAS'), EXP.sample.mosaic = p.ETAS; end
  if isfield(p,'SM'),   EXP.mono.dir = p.SM; end
  if isfield(p,'SS'),   EXP.sample.dir = p.SS; end
  if isfield(p,'SA'),   EXP.ana.dir = p.SA; end
  if isfield(p,'KFIX'), 
    V2K  = 1.58825361e-3;
    VS2E = 5.22703725e-6;
    EXP.Kfixed = p.KFIX; EXP.Lfixed=2*pi/EXP.Kfixed; V=EXP.Kfixed/V2K; EXP.efixed=V*V*VS2E;
  end
  if isfield(p,'FX')
    if p.FX==2, EXP.infin=-1; else EXP.infin=1; end
  end
  if isfield(p,'ALF1') && isfield(p,'ALF2') && isfield(p,'ALF3') && isfield(p,'ALF4')
    EXP.hcol = [ p.ALF1 p.ALF2 p.ALF3 p.ALF4 ];
  end
  if isfield(p,'BET1') && isfield(p,'BET2') && isfield(p,'BET3') && isfield(p,'BET4')
    EXP.vcol = [ p.BET1 p.BET2 p.BET3 p.BET4 ];
  end
  if isfield(p,'AS'), EXP.sample.a=p.AS; end
  if isfield(p,'BS'), EXP.sample.b=p.BS; end
  if isfield(p,'CS'), EXP.sample.c=p.CS; end
  if isfield(p,'AA'), EXP.sample.alpha=p.AA; end
  if isfield(p,'BB'), EXP.sample.beta= p.BB; end
  if isfield(p,'CC'), EXP.sample.gamma=p.CC; end
  if isfield(p,'AX') && isfield(p,'AY') && isfield(p,'AZ')
    EXP.orient1=[p.AX p.AY p.AZ];
  end
  if isfield(p,'BX') && isfield(p,'BY') && isfield(p,'BZ')
    EXP.orient2=[p.BX p.BY p.BZ];
  end
  if isfield(p,'QH'), EXP.QH=p.QH; end
  if isfield(p,'QK'), EXP.QK=p.QK; end
  if isfield(p,'QL'), EXP.QL=p.QL; end
  if isfield(p,'EN'), EXP.W =p.EN; end
  if isfield(p,'H'),  EXP.QH=p.H; end
  if isfield(p,'K'),  EXP.QK=p.K; end
  if isfield(p,'L'),  EXP.QL=p.L; end
  if isfield(p,'W'),  EXP.W =p.W; end
  
% Popovici parameters (pinst: 27)
  if isfield(p,'WB'), EXP.beam.width     =p.WB; end
  if isfield(p,'HB'), EXP.beam.height    =p.HB; end
  if isfield(p,'WS'), EXP.sample.width   =p.WS; end
  if isfield(p,'HS'), EXP.sample.height  =p.HS; end
  if isfield(p,'TS'), EXP.sample.depth   =p.TS; end
  if isfield(p,'WD'), EXP.detector.width =p.WD; end
  if isfield(p,'HD'), EXP.detector.height=p.HD; end
  if isfield(p,'WM'), EXP.mono.width     =p.WM; end
  if isfield(p,'HM'), EXP.mono.height    =p.HM; end
  if isfield(p,'TM'), EXP.mono.depth     =p.TM; end
  if isfield(p,'WA'), EXP.ana.width     =p.WA; end
  if isfield(p,'HA'), EXP.ana.height    =p.HA; end
  if isfield(p,'TA'), EXP.ana.depth     =p.TA; end
  if isfield(p,'L1') && isfield(p,'L2') && isfield(p,'L3') && isfield(p,'L4')
    EXP.arms=[p.L1 p.L2 p.L3 p.L4];
  end
  
  % radius of curvature [m] -> [cm]
  if isfield(p,'RMH'), EXP.mono.rh=100*p.RMH; end
  if isfield(p,'RMV'), EXP.mono.rv=100*p.RMV; end
  if isfield(p,'RAH'), EXP.ana.rh=100*p.RAH; end
  if isfield(p,'RAV'), EXP.ana.rv=100*p.RAV; end
  
  % curvatures from ResCal are in [m-1] -> [cm]
  if isfield(p,'ROMH'), EXP.mono.rh=100/p.ROMH; end
  if isfield(p,'ROMV'), EXP.mono.rv=100/p.ROMV; end
  if isfield(p,'ROAH'), EXP.ana.rh=100/p.ROAH; end
  if isfield(p,'ROAV'), EXP.ana.rv=100/p.ROAV; end

% end ResLibCal_RescalPar2EXP

