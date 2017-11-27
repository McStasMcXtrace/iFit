%   READ_INX (read the INX files -ILL standards)
%
%   Usage: s = read_inx(filename,{instrument})
%                    
%
%         Option: 
%                    instrument = 'in6', 'in5', 'in16', ..., default 'IN5'
%                    to make the difference between the 4-columns and 3-columns
%                    data format.
%                    filename   = filename (inc. filepath if needed)
%
%        Output:
%                 s.header     = All the (string) headers
%                 s.Par        = Double real matrix of Param: Angle, Ei, etc.
%                 s.Mat        = Double real matrix of results (3D).
%
% (c) J. Ollivier, ILL, EUPL 1.1
% See also: read_endf, read_mccode, read_llb_tas, read_anytext


% format spec:
% line 1: Nb tot channels, -, -, -, -, nb channels actually used
% line 2: 'title'
% line 3: Angle, incident energy, transfered wave-vector (Q), mean temperature, -, -
% line 4: -, channel width (musec), -
% [ Energy, Intensity (a.u.), intensity errors ]

%--------------------------------------------------------------------------
%
%  Written by:  JO ~ 2003 for use with plotdata
%  Changes:     JO 2004-2005 for use as a standalone routine
%  Version:     $Date$
%--------------------------------------------------------------------------
function s = read_inx(filename,varargin)


s.header = '';
s.Par    = [];
s.Mat    = [];
%--------------find-path-and-load---------------------------------
if nargin == 0, return; end

FID = fopen(filename);

if FID == -1, 
  error([mfilename ': ERROR: File: ',filename,' not found!']);
end
%----------store-parameters-and-matrix-in-variables----------------
kk=1;
while (1)
   if feof(FID), break;end,
   Dummy      = fgetl(FID); % line 1
   header     = Dummy;
   Dummy      = str2num(Dummy);
   if isempty(Dummy)
     s = [];
     return
   end
   Npoints    = Dummy(length(Dummy));
   Dummy      = fgetl(FID); % line 2= title
   header     = strvcat(header,Dummy);
   Dummy      = fgetl(FID); % line 3
   header     = strvcat(header,Dummy);
   Dummy      = str2num(Dummy);
   Dummy2     = fgetl(FID); % line 4
   header     = strvcat(header,Dummy2);
   s.header(:,:,kk)  =  header; 
   s.Par(kk,:)    = [Dummy,str2num(Dummy2)];
%%   fprintf('Read spectra: %d ...\n',kk); replaced by the "waitbar"
   Old_Cursor = ftell(FID);
   %-----Evaluate-the-number-of-columns-in-file--------
   Dummy      = str2num(fgetl(FID));
   New_Cursor = ftell(FID);
   Ncol       = length(Dummy);
   Ncol       = Ncol(1);
   %-----Rewind-one-step------------------------------
   status     = fseek(FID,Old_Cursor-New_Cursor,'cof');
   %-----Read-one-spectrum----------------------------
   fprintf(1,'%s: file %s, block %i: reading %i x %i values\n', mfilename, filename, kk, Ncol, Npoints);
   Dummy      = fscanf(FID,'%g',Ncol*Npoints);
   for l=1:Npoints,
       dummyM(l,:,kk)=Dummy(1+(l-1)*Ncol:l*Ncol);
   end,
   Dummy = fgetl(FID);
   kk=kk+1;
end,
fclose(FID);
if nargin < 3
   instrument = 'in5';
else
   instrument = varargin{1};
end,   
if strcmp(instrument,'in16')
    if Ncol >= 4
       s.Mat = dummyM(:,2:Ncol,:);
    else
       s.Mat = dummyM;  
    end,
else
    s.Mat = dummyM;
end


