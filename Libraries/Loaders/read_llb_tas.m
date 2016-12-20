function data = read_llb_tas( filename )
%READ_LLB_TAS Reads an LLB/TAS Data set, Computation(model) or Scan file
%   Opens an LLB/TAS data file and returns a data set structure
%   data = read_llb_tas( filename )
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_anytext, read_inx

% built from B. Nennion wfal/read_fich.f
data = [];
if nargin == 0, return; end

[p,f,e]= fileparts(filename);
if isempty(f)
    % this is a dot file (hidden)
    return
end
% the type of LLB TAS data depends on the first letter of the file name, + number
if ~any(upper(f(1)) == 'RCSDP') || ~isfinite(str2double(f(2:end)))
  % error([ mfilename ': ' filename ': ERROR: the filename 1st letter must be any of RCSDP.' ])
  return;
end

if ~isempty(e)
  % the file name must not have extension
  return
end

fid = fopen(filename);
if fid == -1, return; end

% LLB binary format:
% each 'line' is 92 bytes long
recl = 92;

% The file uses a SHORT storage, except for the title
date0      = fread(fid, 5, 'short')'; % the Date [D,M,Y,H,M]
data.date  = datestr([date0([3 2 1 4 5]) 0]);
data.title = deblank(fread(fid, 80, '*char')');

if any(~isstrprop(data.title,'print'))
  fclose(fid);
  disp([ mfilename ': ' filename ': ERROR: the title format is not valid. Probably not an LLB/TAS file format.' ])
  return
end

% the first letter of the filename indicates its type:
switch(upper(f(1)))
    case {'R','C'}
        if upper(f(1)) == 'R'
          data.type = 'Data set';
        else
          data.type = 'Computation/model';
        end
        fseek(fid, (6-1)*recl, 'bof'); % 6-th line of record
        
        % header: [dummy(4),dx(4),xnp,xmon,xif,xk,ctab,nmon,imax,nana,nsp(10)]
        data.hkle  = fread(fid, 4, 'float32');
        data.dhkle = fread(fid, 4, 'float32');
        data.np    = fread(fid, 1, 'float32'); % nb of points in scan
        data.xmon  = fread(fid, 1, 'float32');
        data.fx    = fread(fid, 1, 'float32');
        data.kfix  = fread(fid, 1, 'float32');
        data.nmon  = fread(fid, 1, 'int16');
        data.imax  = fread(fid, 1, 'int16');
        data.idpart= fread(fid, 1, 'int16');
        data.nsp   = fread(fid, 10,'int16');

        % a few tests
        if ~any(data.fx == [1 2]) || data.np < 0 || data.xmon < 0 || data.kfix < 0 ...
          || data.nmon < 0 || data.imax < 0 || data.idpart < 0
          fclose(fid);
          error([ mfilename ': ' filename ': ERROR: the header format is not valid. Probably not an LLB/TAS file format.' ])
        end

        if upper(f(1)) == 'C'
          if data.idpart > 12, data.idpart=2; end
          data.idpart = data.idpart-1;
          % get title or computation index 'no' depending on the 'afit'/'hfit'
          % version
          fseek(fid, (7+data.np-1)*recl, 'bof');
          data.index    = fread(fid, 1, 'int16');
          if ~(1 <= data.index && data.index < 2001)
              % the line after the points is a sub-title: we read it again
              fseek(fid, (7+data.np-1)*recl, 'bof');
              data.subtitle = deblank(fread(fid, 80, '*char')');
          end
        else % R files
          fseek(fid, (7-1)*recl, 'bof'); % 7-th line of record
          data.ns     = fread(fid, 1, 'int16');
          data.angles = fread(fid, 20, 'float32');
          % nang = data.nb;
          % data.angles = fread(fid, 2*data.nb, 'float32');
          % data.angles = data.angles(2:2:end);
        end
        
        % each data line:
        kmax = 12+max(2,data.imax)+data.idpart;
        % imm=12+imax;
        % [short=no][real(1:kmax), (short)nsp1, (real)ture(1:nsp1)]
        % [short=no][(real)dx][real(1:(kmax-5))][real=cnt][real=y]
        data.signal = zeros(kmax, data.np);
        for i=1:data.np
          fseek(fid, (6+i-1)*recl, 'bof');
          data.index(i)=fread(fid, 1, 'int16');
          data.signal(:,i) = fread(fid, kmax, 'float32');
        end
        data.signal= data.signal';
        data.index = data.index';

    case 'S'
        data.type = 'Scan';
        % spectrometer ID
        fseek(fid, (2-1)*recl, 'bof'); % 2-nd line of record
        data.spec = fread(fid, 12, 'int16');

        fseek(fid, (3-1)*recl, 'bof'); % 3-rd line of record
        data.eps  = fread(fid, 9, 'int16');

        nspec = 1; % should get files for each spectrometer to know how to get
                   % the proper ID
                   
        % configurations
        switch data.spec
        case 1
          % 1T/2T
          vars = {'M1 ','M2 ','E1 ','E2 ','A1 ','A2 ','A7 ', ...
          'A8 ','GM ','MON','GI ','GS ','FOC','CA ','CM2','TM ','OHM',' K', ...
          ' H ','','','','','','GM ','TM ','CM1','CM2','MON','','','','','','','','d1r','d1l', ...
          'd1d','d1u','d2r','d2l','d2d','d2u','','','','', ...
          'M1 ','M2 ','M3 ','M4 ','E1 ','E2 ','A1 ','A2 '};
        case 2
          % 4F1
          vars = {'GM1','GM2','GI ','GS ','FOC','CA ','B1 ','B2 ','OHM',' K',' H ', ...
           '','','','','','','','','','','','','','','','','', ...
           'd1r','d1l','d1d','d1u','d2r','d2l','d2d','d2u','rco', 'vco','','',...
           'M1 ','M2 ','M3 ','M4 ','E1 ','E2 ','A1 ','A2 '};
        case 3
          % 4F2
          vars = {'GM1','GM2','GI ','GS ','FOC','CA ',' M ',' C ','OHM',' K',' H ', ...
           '','','','','','','','','','','','','','','','','', ...
           'd1r','d1l','d1d','d1u','d2r','d2l','d2d','d2u','rco', 'vco','','', ...
           'M1 ','M2 ','M3 ','E1 ','E2 ','A1 ','A2 ','D1 '};
        case 4
          % G43
          vars = {'D2 ','D3 ','GI ','GS ','FOC','CA ','A1 ','A2 ','OHM','K',' H '};
        end

        % header: [nb angles, np]
        fseek(fid, (6-1)*recl, 'bof'); % 6-th line of record
        data.nb   = fread(fid, 1, 'int16'); data.np  = fread(fid, 1, 'int16');
        data.nmon = fread(fid, 1, 'int16'); data.mon = fread(fid, 1, 'int16');
        data.nana = fread(fid, 1, 'int16'); data.nsp = fread(fid, 10, 'int16');

        

        fseek(fid, (7-1)*recl, 'bof'); % 7-th line of record
        data.ns  = fread(fid, 1, 'int16');
        % nang = data.nb;
        data.angles = fread(fid, 2*data.nb, 'float32');
        data.angles = data.angles(2:2:end);
        
        for i=1:data.np
          fseek(fid, (6+i-1)*recl, 'bof');
          % each line is [no,id,angle] (nb-1)*[id,cont] [xmon,ymon,cnt]
          data.index(i)= fread(fid, 1, 'int16');
          tmp       = fread(fid, 2*data.nb, 'float32');
          angles    = tmp(2:2:end);
          id        = tmp(1:2:end);
          tmp       = fread(fid, 3, 'float32');
          xt = tmp(1); xm=tmp(2); xx = tmp(3);
          nsp       = fread(fid, 1, 'int16');
          ture      = fread(fid, nsp, 'float32');
          % the scanned parameters are: vars{1:nb} with values: angles(1:nb)
          % the temperatures are:       ture(1:nsp)
          % counting is:                xx
          for j=1:data.nb
            name = strtrim(vars{j});
            if isempty(name), continue; end
            data.(name)(i) = angles(j); % update vector
          end
          data.signal(i) = xx;
          for j=1:nsp
            name = [ 'T_' num2str(j) ];
            data.(name)(i) = ture(j);
          end
          
        end
        data.columns = vars;
    case 'D'
        data.type = 'Counting';
    otherwise
        data = [];
end

% test the date [D,M,Y,H,M]
D = date0(1); M=date0(2); Y=date0(3); H=date0(4); m=date0(5);
if any([ D M Y H m ] <= 0) || D > 31 || M > 12 || m > 60
  data.date = date;
  disp([ 'WARNING: ' mfilename ': ' f e ': the date format is not valid. Using today.' ])
end

% close the file
fclose(fid);

end

