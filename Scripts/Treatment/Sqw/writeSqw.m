function [this, parameters] = writeSqw(this, parameters)
% [this, parameters] = writeSqw(this, parameters) : write Sqw file for McStas 
%   from NetCDF nMoldyn output
%
%   Metadata may be given in the parameters input argument, and then appears in 
%   the Sqw file header (physical parameters, see below).
%   This function requires iFit to be installed [ifit.mccode.org]
%
% Input:
%   this:       NetCDF nMoldyn SQW filename (DCSF, DISF) or iData object
%   parameters: vector or structure containing the simulation parameters
%   parameters.density     [g/cm3] Material density
%   parameters.weight      [g/mol] Material weight
%   parameters.T_m         [K] Melting temperature
%   parameters.T_b         [K] Boiling temperature
%   parameters.MD_at       Number of atoms in the molecular dynamics simulation
%   parameters.MD_box      [Angs] Simulation box size
%   parameters.Temperature [K] Temperature
%   parameters.dT          [K] Temperature accuracy
%   parameters.MD_duration [ps] Molecular dynamics duration
%   parameters.D           [cm^2/s] Diffusion coefficient
%   parameters.sigma_coh   [barn] Coherent scattering neutron cross section
%   parameters.sigma_inc   [barn] Incoherent scattering neutron cross section
%   parameters.sigma_abs   [barn] Absorption neutron cross section
%   parameters.c_sound     [m/s] Sound velocity. 'v_sound' is also supported.
%   parameters.At_number   Atomic number Z
%   parameters.Pressure    [bar] Material pressure
%   parameters.Material    Material description
%   parameters.Phase       Material state
%   parameters.Scattering  scattering process
%   parameters.Lambda      [Angs] Neutron wavelength used for the measurement
%   parameters.Instrument neutron spectrometer used for measurement
%   parameters.Comment     any additional comment/description
%   parameters.C_s         [J/mol/K]   heat capacity
%   parameters.Sigma       [N/m] surface tension
%   parameters.Eta         [Pa s] viscosity
%   parameters.L           [J/mol] latent heat
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
if nargin < 2, parameters=[]; end

if isempty(this), return; end

filename = [ this.Source '.sqw' ];

fields={'density [g/cm3] Material density'	'weight [g/mol] Material weight'	'T_m [K] Melting T'	...
    'T_b [K] Boiling T'	'MD_at Number of atoms in the molecular dynamics simulation'	...
    'MD_box [Angs] Simulation box size'	'Temperature [K]'	'dT [K] T accuracy'...
    'MD_duration [ps] Molecular dynamics duration'	'D [cm^2/s] Diffusion coefficient'...
    'sigma_coh [barn] Coherent scattering neutron cross section'...
    'sigma_inc [barn] Incoherent scattering neutron cross section'	...
    'sigma_abs [barn] Absorption neutron cross section'	'c_sound [m/s] Sound velocity'...
    'At_number Atomic number Z'};

if ~isempty(parameters) && isnumeric(parameters)
  if length(parameters) == length(fields)+1
    fields{end}='multiplicity  [atoms/unit cell] number of atoms/molecules per scattering unit cell';
  elseif length(parameters) ~= length(fields)
    error(sprintf('%s: Error: %s: length of parameters=%i length of fields=%i', ...
      mfilename, filename, length(parameters), length(fields)));
  end
  p=parameters; parameters=[];
  
  for index=1:length(fields)
    name = strtok(fields{index});
    parameters.(name) = p(index);
  end
else
  % build a parameter structure by searching field names from 'fields' in the object
  fields{end+1} = 'Pressure [bar] Material pressure';
  fields{end+1} = 'v_sound [m/s] Sound velocity';
  fields{end+1} = 'Material';
  fields{end+1} = 'Phase Material state';
  fields{end+1} = 'Scattering process';
  fields{end+1} = 'Lambda [Angs] Neutron wavelength used for the measurement';
  fields{end+1} = 'Wavelength [Angs] measurement neutron wavelength';
  fields{end+1} = 'Instrument neutron spectrometer used for measurement';
  fields{end+1} = 'Comment';
  fields{end+1} = 'density [g/cm3] Material density ';
  fields{end+1} = 'C_s [J/mol/K]   heat capacity';
  fields{end+1} = 'Sigma [N/m] surface tension';
  fields{end+1} = 'Eta [Pa s] viscosity';
  fields{end+1} = 'L [J/mol] latent heat';
  fields{end+1} = 'classical [0=from measurement, with Bose factor included, 1=from MD, symmetric]';
  for index=1:length(fields)
    name = strtok(fields{index});
    if ~isfield(parameters, name)
      if isfield(this, name)
        parameters.(name) = get(this, name);
      elseif isfield(this.Data, name)
        parameters.(name) = this.Data.(name);
      elseif isfield(this.Data, name)
        parameters.(name) = this.Data.(name);
      else
        parameters.(name) = NaN;
      end
    end
  end
end

% add parameters from the jobinfo, when applicable
if isfield(this, 'title_nc')
  t = strtrim(this.title_nc); t(~isstrprop(t,'print')) = ' ';
  parameters.MD_type = strtrim(t);
  fields{end+1} = 'MD_type Type of nMoldyn processing';
  parameters.Scattering = strtok(parameters.MD_type,' _');
end
if isfield(this, 'history_nc')
  [dummy, parameters.Trajectory] = strtrim(strtok(this.history_nc));
  [p,f,e] =fileparts(parameters.Trajectory);
  parameters.Trajectory = [ f e ];
  fields{end+1} = 'Trajectory Molecular dynamics trajectory';
end

if isfield(this, 'jobinfo')
  if isfield(this.jobinfo, 'frequencyunits')
    if ~strcmp(this.jobinfo.frequencyunits,'meV')
      error(sprintf(1,'%s: Warning: %s: frequency is not given in [meV]', mfilename, filename));
    end
  end
  if isfield(this.jobinfo, 'qunits')
    if ~strcmp(this.jobinfo.qunits,'ang^-1')
      error(sprintf(1,'%s: Warning: %s: momentum is not given in [ang^-1]', mfilename, filename));
    end
  end
  if isfield(this.jobinfo, 'resolution')
    parameters.MD_resolution = this.jobinfo.resolution;
    fields{end+1} = 'MD_resolution [meV] Energy convolution width used in the FFT (r,t)->(q,w)';
  end
  if isfield(this.jobinfo, 'weight')
    parameters.Scattering = this.jobinfo.weights;
  end
  if isfield(this.jobinfo, 'trajectory') && ~isfield(parameters, 'Trajectory')
    parameters.Trajectory = this.jobinfo.trajectory;
    [p,f,e] =fileparts(parameters.Trajectory);
    parameters.Trajectory = [ f e ];
    fields{end+1} = 'Trajectory Molecular dynamics trajectory';
  end
  if isfield(this.jobinfo, 'output')
    parameters.MD_sqw = fliplr(strtok(fliplr(strtrim(this.jobinfo.output)), ': /\'));
    [p,f,e] =fileparts(parameters.MD_sqw);
    parameters.MD_sqw = [ f e ];
    fields{end+1} = 'MD_sqw Molecular dynamics dynamic structure factor file from nMoldyn';
  end
end

if isfield(parameters, 'Scattering')
  if strfind(lower(parameters.Scattering),'incoherent')
    parameters.Scattering = 'incoherent';
  elseif strfind(lower(parameters.Scattering),'coherent')
    parameters.Scattering = 'coherent';
  else
    parameters.Scattering = 'total';
  end
end

if isfield(parameters, 'Trajectory')
  parameters.Material = fliplr(strtok(fliplr(strtrim(parameters.Trajectory)), ': /\'));
  [dummy,parameters.Material] = fileparts(parameters.Material);
  if parameters.Material(1) == 'l'
    parameters.Phase = 'liquid';
  elseif parameters.Material(1) == 'p'
    parameters.Phase = 'powder';
  end
  if parameters.Material(1) == 'l' || parameters.Material(1) == 'p'
    parameters.Material=parameters.Material(2:length(parameters.Material)); 
  end
else
  parameters.Trajectory = NaN;
end

% check/set Material, Phase, Scattering and Pressure
if all(isnan(parameters.Phase))
  if ~isempty(findstr(this,'liquid')) || ~isempty(strfind(this.Source,'liq'))
    parameters.Phase = 'liquid';
  elseif ~isempty(findstr(this,'powder')) || ~isempty(strfind(this.Source,'pow'))
    parameters.Phase = 'powder';
  elseif ~isempty(findstr(this,'gas')) || ~isempty(strfind(this.Source,'gas'))
    parameters.Phase = 'gas';
  elseif ~isnan(parameters.Pressure)
    parameters.Phase = [ num2str(parameters.Pressure) 'bar' ];
  else
    parameters.Phase = 'phase';
  end
end

if all(isnan(parameters.Material))
  parameters.Material = 'material';
end

this.parameters=parameters;

% check orientation of Signal
sqw=this{0}; sqw(isnan(sqw))=0;
w=this{2}; w=w(:)'; 
q=this{1}; q=q(:)'; 
if size(sqw,1) == length(w) && size(sqw,2) == length(q)
  sqw=sqw';
end

% transpose the object axis and Signal
%iw=find(abs(w) < 150); w=w(iw);
%sqw=sqw(:,iw);

% format filename
[p,f] = fileparts(filename);
% p = '/home/farhi/Simulations/Sqw-lib/Sqw';
p = pwd;
if ~isdir(p), p=pwd; end
f = sprintf('%s_%s_%i_%s.sqw', ...
  parameters.Material, parameters.Phase(1:3), ceil(parameters.Temperature), parameters.Scattering(1:3));
filename = fullfile(p,f);


% write the file header
  [fid,message] = fopen(filename, 'w+');
  if fid < 0
    error([ mfilename ': can not open file ' filename ' because ' message ]); 
  else
    fprintf(1,'Opening %s\n', filename);
  end

  fprintf(fid,'# Format: Sqw data file for Isotropic_Sqw <http://www.mcstas.org>\n');
  
  fprintf(fid,'# %s %s', parameters.Phase, parameters.Material);
  if all(~isnan(parameters.MD_at)) && all(~isnan(parameters.MD_duration)) && all(~isnan(parameters.MD_box))
    fprintf(fid,' molecular dynamics simulation using %g atoms for %g [ps]. box=%g Angs\n', ...
      parameters.MD_at, parameters.MD_duration, parameters.MD_box);
  else
    fprintf(fid, '\n');
  end
    
  fprintf(fid,'# Sqw data from E. Farhi');
  if all(~isnan(parameters.Trajectory))
    fprintf(fid, ', obtained from nMoldyn <http://forge.ill.fr/projects/nmoldyn>\n');
  else
    fprintf(fid, '\n');
  end
  fprintf(fid,'#\n');
  
  fprintf(fid,'# title: %s %s: S(q,w) %s part\n', parameters.Phase, parameters.Material, parameters.Scattering);

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

% write the q,w,Sqw data
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

parameters.filename = filename;
this.parameters=parameters;

