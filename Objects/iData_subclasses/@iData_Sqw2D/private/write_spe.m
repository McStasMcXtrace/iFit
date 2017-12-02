function filename = write_spe(spe, filename)

if nargin < 2,        filename = []; end
if isempty(filename), filename = [ 'iFit_' spe.Tag '.spe' ]; end

% INX format:
% s.header     = All the (string) headers
%   line 1: Nb tot channels, -, -, -, -, nb channels actually used
%   line 2: 'title'
%   line 3: Angle, incident energy, transfered wave-vector (Q), mean temperature, -, -
% line 4: -, channel width (musec), -
% s.Par        = Double real matrix of Param: Angle, Ei, etc.
% s.Mat        = Double real matrix of results (3D).

phi = getaxis(spe, 2);
w   = getaxis(spe, 1);
signal = getaxis(spe, 0);
err    = getaxis(spe, 'Error');

[FID,mess] = fopen(filename,'w');
if FID == -1, 
  disp(mess);
  error([mfilename ': ERROR: File ' filename ' can not be opened.']);
end

% header
F = '%4d\t';
fprintf(FID,[F,F,'\n'], size(spe, 2), size(spe, 1));  % [nbangles nbchan]

% angles
fprintf(FID,'### Phi Grid\n');
phi = mean(phi,1);  % mean angles along rows
str = num2str(phi);
str(:,end+1) = sprintf('\n');
str = str';
str = str(:)';
fprintf(FID, '%s', str);

% energies
fprintf(FID,'### Energy Grid\n');
w   = mean(w,2);    % mean energies along columns
str = num2str(w);
str(:,end+1) = sprintf('\n');
str = str';
str = str(:)';
fprintf(FID, '%s', str);

% S(phi,w)
fprintf(FID,'### S(Phi,w)\n');
str = num2str(signal);
str(:,end+1) = sprintf('\n');
str = str';
str = str(:)';
fprintf(FID, '%s', str);

% errors
fprintf(FID,'### Errors\n');
str = num2str(err);
str(:,end+1) = sprintf('\n');
str = str';
str = str(:)';
fprintf(FID, '%s', str);

fclose(FID);
