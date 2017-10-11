function filename = Sqw_McStas(this, filename)
% [this, parameters] = Sqw_McStas(this) : write Sqw file for McStas 
%   from NetCDF nMoldyn output
%
%   Metadata may be given as aliases in the object, and then appears in 
%   the Sqw file header (physical parameters).
%   This function requires iFit to be installed [ifit.mccode.org]
%
% Input:
%   this:       NetCDF nMoldyn SQW filename (DCSF, DISF) or iData object
%   filename:   output file to be written for McStas. Standard naming is:
%               Material_Phase_Temperature_Scattering where Phase and Scattering are
%               3 chars (e.g. liq,pow,gas and coh,inc,tot).
%
% (c) E. Farhi ILL/DS/CS 2012.

if nargin < 1, this=[]; end
if ischar(this),  
  [p,f,e] = fileparts(this);
  if isempty(e), this = [ this '.nc' ]; end
  if ~any(strncmp(this,{'DISF','DCSF'}, 4))
    [this_i,p_i] = feval(mfilename, [ 'DISF_' this ], parameters);
    [this_c,p_c] = feval(mfilename, [ 'DCSF_' this ], parameters);
    this       = [ this_c, this_i ];
    parameters = { p_c,    p_i };
    return
  else
    this=iData(this);
  end
elseif isempty(this), this=iData(''); end

if ~isa(this, 'iData')
  disp([ mfilename ': ERROR: The data set should be an iData object, and not a ' class(this) ]);
  filename = [];
  return; 
end

% check object and orientation of Signal
this = Sqw_check(this);
if isempty(this), return; end

if nargin < 2, filename = ''; end

% call Sqw_parameters so that we have the final parameters and the comments
[this,parameters,fields] = Sqw_parameters(this);

sqw=this{0}; sqw(isnan(sqw))=0;
w=this{1}; w=w(:)'; 
q=this{2}; q=q(:)'; 

% format filename
if isempty(filename)
  if ~isdir(this.Source)
    [~,filename] = fileparts(this.Source);
  else
    filename = [ 'Sqw_' b.Tag ];
  end
  filename = [ filename '.sqw' ];
end

% write the file header --------------------------------------------------------
[fid,message] = fopen(filename, 'w+');
if fid < 0
  error([ mfilename ': can not open file ' filename ' because ' message ]); 
else
  fprintf(1,'Opening %s\n', filename);
end

fprintf(fid,'# Format: Sqw data file for Isotropic_Sqw <http://www.mcstas.org>\n');

if isfield(parameters,'Phase') && isfield(parameters,'Material')
  fprintf(fid,'# %s %s', parameters.Phase, parameters.Material);
else
  fprintf(fid,'#\n');
end
if isfield(parameters,'MD_at') && isfield(parameters,'MD_duration') && isfield(parameters,'MD_box')
  fprintf(fid,' molecular dynamics simulation using %g atoms for %g [ps]. box=%g Angs\n', ...
    parameters.MD_at, parameters.MD_duration, parameters.MD_box);
else
  fprintf(fid, '\n');
end
if isfield(parameters,'Trajectory')
  fprintf(fid, ', obtained from nMoldyn/MDANSE <http://forge.ill.fr/projects/nmoldyn>\n');
else
  fprintf(fid, '\n');
end
fprintf(fid,'#\n');

if isfield(parameters,'Phase') && isfield(parameters,'Material') && isfield(parameters,'Scattering')
  fprintf(fid,'# title: %s %s: S(q,w) %s part\n', ...
    parameters.Phase, parameters.Material, parameters.Scattering);
end

fprintf(fid,'# Date: %s\n', datestr(now));
[p,f,e] =fileparts(this.Source);
fprintf(fid,'# Source: %s%s\n', f,e);
[p,f,e] =fileparts(filename);
fprintf(fid,'# filename: %s%s\n', f,e);
fprintf(fid,'# format: Sqw data file for Isotropic_Sqw (McStas)\n');
i = ~isnan(sqw);
fprintf(fid,'# signal: Min=%g; Max=%g; Mean=%g; sum=%g;\n', min(sqw(i)), max(sqw(i)), mean(sqw(i)), sum(sqw(i)));
fprintf(fid,'# type: array_2d(%i,%i)\n', length(q), length(w));
fprintf(fid,'# xylimits: %g %g %g %g\n', min(q), max(q), min(w), max(w)); 
fprintf(fid,'# xlabel: Wavevector [Angs-1]\n'); 
fprintf(fid,'# ylabel: Energy [meV]\n'); 
fprintf(fid,'#\n');

if isstruct(parameters)
  parameters = orderfields(parameters);
  orderedfields=fieldnames(parameters);
  fprintf(fid,'# Physical parameters:\n');
  for index=1:length(orderedfields)
    field_found=0;
    for field_index=1:length(fields)
      if strcmpi(orderedfields{index}, strtok(fields(field_index)))
        field_found=1;
        break;
      end
    end
    if ~field_found, continue; end
    [name, comment] = strtok(fields{field_index});
    comment = strtrim(comment);
    if isfield(parameters, name)
      if isnumeric(parameters.(name))
        if isfinite(parameters.(name))
          fprintf(fid, '# %-13s %-8.3g %s\n', name, parameters.(name), comment);
          parameters.(name)=NaN;
        end
      else
        t=parameters.(name); t(~isstrprop(t,'print')) = ' '; parameters.(name)=t;
        fprintf(fid, '# %-13s %-s %s\n', name, parameters.(name), comment);
        parameters.(name)=NaN;
      end
    end
  end
end

% write the q,w,Sqw data -------------------------------------------------------

% add q text
fprintf(fid, '#\n');
fprintf(fid, '# WAVEVECTOR vector of m=%i values %g:%g in [Angstroem-1]: q\n',length(q), min(q),max(q));
% add q
fprintf(fid, '%g ', q);

% add w text
fprintf(fid, '\n');
fprintf(fid, '# ENERGY vector of n=%i values %g:%g in [meV]: w\n',length(w),min(w),max(w));
% add w
fprintf(fid, '%g ', w);

% add sqw text
fprintf(fid, '\n');
fprintf(fid, '# matrix of S(q,w) values (m=%i rows x n=%i columns), one row per q value: sqw\n', ...
  length(q),length(w));
% add sqw
for k=1:length(q)
  fprintf(fid, '%s', num2str(sqw(k,:)));
  fprintf(fid, '\n');
end
fprintf(fid, '# end of Sqw file %s\n', filename);
fclose(fid);
fprintf(1,'DONE %s\n', filename);


