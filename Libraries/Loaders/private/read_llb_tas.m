function data = read_llb_tas( filename )
%READ_LLB_TAS Reads an LLB/TAS Data set, Computation(model) or Scan file
%   Opens an LLB/TAS data file and returns a data set structure
data = [];

fid = fopen(filename);
if fid == -1, return; end

% LLB binary format:
% each 'line' is 92 bytes long
recl = 92;

% The file uses a SHORT storage, except for the title
date0      = fread(fid, 5, 'short')'; % the Date [D,M,Y,H,M]
data.date  = datestr([date0([3 2 1 4 5]) 0]);
data.title = deblank(fread(fid, 80, '*char')');

[p,f,e]= fileparts(filename);

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
          error([ mfilename ': ' filename ': ERROR: the header format is not valid. Probably not an LLB/TAS file format.' ])
        end

        if upper(f(1)) == 'C'
          if data.idpart > 12, data.idpart=2; end
          data.idpart = data.idpart-1;
          % get title or computation index 'no' depending on the 'afit'/'hfit'
          % version
          fseek(fid, (7+data.np-1)*recl, 'bof');
          data.no    = fread(fid, 1, 'int16');
          if ~(1 <= data.no && data.no < 2001)
              % the line after the points is a sub-title: we read it again
              fseek(fid, (7+data.np-1)*recl, 'bof');
              data.subtitle = deblank(fread(fid, 80, '*char')');
          end
        end
        
        % each data line:
        kmax = 12+max(2,data.imax)+data.idpart;
        % imm=12+imax;
        % [short=no][real(1:kmax), (short)nsp1, (real)ture(1:nsp1)]
        % [short=no][(real)dx][real(1:(kmax-5))][real=cnt][real=y]
        data.array = zeros(kmax, data.np);
        for i=1:data.np
          fseek(fid, (6+i-1)*recl, 'bof');
          data.no(i)=fread(fid, 1, 'int16');
          data.array(:,i) = fread(fid, kmax, 'float32');
        end
        data.array = data.array';
        data.no    = data.no';

    case 'S'
        data.type = 'Scan';
        % header: [nb, np]
        data.nb = block(1); data.np = block(2); block(1:2)=[];
        % each ligne is [no,id,x] (nb-1)*[id,cont] [xmon,ymon,cnt]
        data.data = reshape(block, data.np, numel(block)/data.np);
    case 'D'
        data.type = 'Counting';
    otherwise
        data = [];
end

% test the date [D,M,Y,H,M]
D = date0(1); M=date0(2); Y=date0(3); H=date0(4); m=date0(5);
if any([ D M Y H m ] <= 0) || D > 31 || M > 12 || m > 60
  data.date = date;
  warning([ mfilename ': ' filename ': WARNING: the date format is not valid. Using today.' ])
end
if any(~isstrprop(data.title,'print'))
  data.title
  error([ mfilename ': ' filename ': ERROR: the title format is not valid. Probably not an LLB/TAS file format.' ])
end
% the type of LLB TAS data depends on the first letter of the file name
if ~any(upper(f(1)) == 'RCSDP')
  error([ mfilename ': ' filename ': ERROR: the filename 1st letter must be any of RCSDP.' ])
  return;
end

% close the file
fclose(fid);

end

