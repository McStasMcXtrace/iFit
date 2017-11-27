function filename = write_inx(spe, filename)
% write_inx: write an INX file from an S(phi, w) data set

if nargin < 2,        filename = []; end
if isempty(filename), filename = [ 'iFit_' spe.Tag '.inx' ]; end

% INX format:
% s.header     = All the (string) headers
%   line 1: Nb tot channels, -, -, -, -, nb channels actually used
%   line 2: 'title'
%   line 3: Angle, incident energy, transfered wave-vector (Q), mean temperature, -, -
%   line 4: -, channel width (musec), -
% Example lines 1-4:
%   size(spe,1)+3    1    2    0    0    0    0 size(spe,1)
%   Wood  MP  RNase A wet 305K DOS         
%   Angle(index)   Energy  2*pi/Wavelength  Temperature   0.0  -2.0
%                  0.0000  chan_width_t  0.0000
% s.Par        = Double real matrix of Param: Angle, Ei, etc.
% s.Mat        = Double real matrix of results (3D).

% we have a 2D input [Energy,Angle]
spe = meshgrid(spe);    % make it a square grid, spe=Sqw_q2phi(s)
w   = getaxis(spe, 1);
phi = getaxis(spe, 2);
sig = getaxis(spe, 0);
err = getaxis(spe, 'Error');
titl= title(spe);

% get the parameters we need from 'parameters' and compute the others
[spe,lambda,distance,chwidth,energy,wavevector] = Sqw_search_lambda(spe);

if isempty(chwidth), chwidth = 1.0; end
if isempty(distance),distance= 4.0; end
Temperature=Sqw_getT(spe);        % [K]         Temperature

Npoints = size(spe,1);

[FID,mess] = fopen(filename,'w');
if FID == -1, 
  disp(mess);
  error([mfilename ': ERROR: File ' filename ' can not be opened.']);
end

for index=1:size(spe, 2)        % loop on nb of blocks/angle channels
  if isvector(phi), Angle = phi(index)
  else              Angle = mean(phi(:,index)); end
  % line 1
  fprintf(FID, '%s\n', num2str([Npoints+3 1 2 0 0 0 0 Npoints]));
  
  % line 2
  fprintf(FID, '%s\n', titl);
  
  % line 3
  fprintf(FID, '%s\n', num2str([ Angle energy wavevector Temperature 0.0 -2.0], '%.4e'));
  
  % line 4
  fprintf(FID, '               0.0000  %.4e  0.0000\n', chwidth);
  
  % matrix [ Energy, Intensity (a.u.), intensity errors ]
  matrix = [ w(:,index) sig(:,index) err(:,index) ];
  str = num2str(matrix);
  str(:,end+1) = sprintf('\n');
  str = str';
  str = str(:)';
  fprintf(FID, '%s', str);
end
fclose(FID);
