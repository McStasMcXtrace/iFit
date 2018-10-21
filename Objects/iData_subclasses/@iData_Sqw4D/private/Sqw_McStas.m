function filename = Sqw_McStas(this, filename)
% [this, parameters] = Sqw_McStas(this) : write Sqw 4D file for McStas 
%
%   Metadata may be given as aliases in the object, and then appears in 
%   the Sqw file header (physical parameters).
%   This function requires iFit to be installed [ifit.mccode.org]
%
% Input:
%   this:       SQW iData object
%   filename:   output file to be written for McStas. Standard naming is:
%               Material_Phase_Temperature_Scattering where Phase and Scattering are
%               3 chars (e.g. liq,pow,gas and coh,inc,tot).
%
% (c) E. Farhi ILL/DS/CS 2012.


if ~isa(this, 'iData')
  disp([ mfilename ': ERROR: The data set should be an iData object, and not a ' class(this) ]);
  filename = [];
  return; 
end

% check object and orientation of Signal
% this = Sqw_check(this); % done in iData_Sqw2D
if isempty(this), return; end
if nargin < 2, filename = ''; end

% we must generate a [ h k l w S ] event list
this = meshgrid(this, 'grid');
this = event(this);
h=getaxis(this,1); 
k=getaxis(this,2);
l=getaxis(this,3);
w=getaxis(this,4);
sqw=getaxis(this,0);

% call Sqw_parameters so that we have the final parameters and the comments
parameters = parseparams(this);

% format filename
if isempty(filename)
  if ~isdir(this.Source)
    [~,filename] = fileparts(this.Source);
  else
    filename = [ 'Sqw4_' this.Tag ];
  end
  filename = [ filename '.sqw4' ];
end
[~,~,e] = fileparts(filename);
if isempty(e), filename = [ filename '.sqw4' ]; end

% write the file header --------------------------------------------------------
[fid,message] = fopen(filename, 'w+');
if fid < 0
  error([ mfilename ': can not open file ' filename ' because ' message ]); 
else
  fprintf(1,'Opening %s\n', filename);
end

fprintf(fid,'# Format: Sqw data file (4D) for Single_crystal_inelastic <http://www.mcstas.org>\n');
fprintf(fid,'# title: %s\n', this.Title);
parameters.Phase = 'Single crystal';
if isfield(parameters,'Phase') && isfield(parameters,'Material')
  fprintf(fid,'# %s %s\n', parameters.Phase, parameters.Material);
end

if isfield(parameters,'Phase') && isfield(parameters,'Material') && isfield(parameters,'Scattering')
  fprintf(fid,'# %s %s: S(q,w) %s part\n', ...
    parameters.Phase, parameters.Material, parameters.Scattering);
end
% check for lattice parameters
if isfield(parameters, 'cell') && numel(parameters.cell) == 6
  if ~isfield(parameters, 'lattice_a')  parameters.lattice_a = cell(1); end
  if ~isfield(parameters, 'lattice_b')  parameters.lattice_b = cell(2); end
  if ~isfield(parameters, 'lattice_c')  parameters.lattice_c = cell(3); end
  if ~isfield(parameters, 'lattice_aa') parameters.lattice_aa = cell(4); end
  if ~isfield(parameters, 'lattice_bb') parameters.lattice_bb = cell(5); end
  if ~isfield(parameters, 'lattice_cc') parameters.lattice_cc = cell(6); end
end

fprintf(fid,'#\n');
fprintf(fid,'# Date: %s\n', datestr(now));
[p,f,e] =fileparts(this.Source);
fprintf(fid,'# Source: %s%s\n', f,e);
[p,f,e] =fileparts(filename);
fprintf(fid,'# filename: %s%s\n', f,e);
fprintf(fid,'# format: Sqw data file (4D) for Single_crystal_inelastic (McStas)\n');
fprintf(fid,'# signal: Min=%g; Max=%g; Mean=%g; sum=%g;\n', min(sqw(:)), max(sqw(:)), mean(sqw(:)), sum(sqw(:)));
fprintf(fid,'# type: event list(%i)\n', length(sqw));
fprintf(fid,'#\n');

if isstruct(parameters)
  parameters = orderfields(parameters);
  orderedfields=fieldnames(parameters);
  fprintf(fid,'# Physical parameters:\n');
  for index=1:length(orderedfields)

    name    = orderedfields{index};
    comment = label(this, name);
    %[name, comment] = strtok(fields{field_index});
    comment = strtrim(comment);
    t=comment; t(~isstrprop(t,'print')) = ' '; comment=t;
    if isfield(parameters, name)
      try
        val = get(this, name);
      catch
        val = parameters.(name);
      end
      if isnumeric(val) && numel(val) > 1
        val = mat2str(val);
      end
      if isnumeric(val)
        if isfinite(val)
          fprintf(fid, '# %-13s %-8.3g %s\n', name, val, comment);
          parameters.(name)=NaN;
        end
      else
        if     iscell(val) && ischar(val{1}), val = sprintf('%s ', val{:});
        elseif iscell(val) && isnumeric(val{1}), val = sprintf('%f ', val{:});
        end
        t=val; t(~isstrprop(t,'print')) = ' '; val=t;
        fprintf(fid, '# %-13s %-s %s\n', name, val, comment);
        parameters.(name)=NaN;
      end
    end
  end
end

% write the h,k,l,w,Sqw data -------------------------------------------------------

% add text
fprintf(fid, '#\n');
fprintf(fid, '# File Format:\n');
fprintf(fid, '# column_h 1\n');
fprintf(fid, '# column_k 2\n');
fprintf(fid, '# column_l 3\n');
fprintf(fid, '# column_E 4\n');
fprintf(fid, '# column_S 5\n');
fprintf(fid, '#\n');
fprintf(fid, '# limits: h=[%g:%g] k=[%g:%g] l=[%g:%g] w=[%g:%g]\n', ...
  min(h), max(h), min(k), max(k), min(l), max(l), min(w), max(w));
fprintf(fid, '# h k l En S(q,w)\n');
str = num2str([ h k l w sqw ], '%g ');
str(:,end+1) = sprintf('\n');
str = str';
str = str(:)';
fprintf(fid, '%s', str);
fprintf(fid, '# end of Sqw4 file %s\n', filename);
fclose(fid);
fprintf(1,'DONE %s\n', filename);

