function [s, header, data] = read_hbin(filename)
% read ILL Cyclops image
%   s = read_hbin(filename)
%
% The HBIN format is an extended TIFF format.
%  <https://www.ill.eu/instruments-support/instruments-groups/instruments/cyclops/description/instrument-layout/>
%
% (c) E.Farhi, ILL. License: EUPL.
% See also: read_edf, read_adsc, read_cbf, read_sif, read_mar, read_spe, read_fits, read_image

s = [];
if nargin == 0, return; end

fid = fopen(filename);
if fid == -1, return; end

% read length of header
offset = fgetl(fid);
offset = str2struct(offset);

if isfield(offset,'Image_Offset')
  offset = str2double(strtok(offset.Image_Offset));
else offset = [];
end

% read header length
header = fgetl(fid);
header = str2struct(header);

if isfield(header,'Header_Length')
  header = str2double(strtok(header.Header_Length));
else header = [];
end

if ~isempty(offset) && isfinite(offset) && ~isempty(header) && isfinite(header)

  % reset file
  fseek(fid, 0, 'bof');

  header = fread(fid, header, 'uint8=>char');  % the header
  % convert header to structure
  s = str2struct(header');
  
  fseek(fid, offset, 'bof');
  s.Signal   = uint16(fread(fid, Inf, 'uint16'));       % to end
  
  if isfield(s, 'Nrow_Ncol'), s.Signal = reshape(s.Signal, s.Nrow_Ncol'); end
end
fclose(fid);

