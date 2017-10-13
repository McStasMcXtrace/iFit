function data = read_opus(filename)
% reads an Opus Brucker file entirely.
%  data = read_opus(filename)
%
% all sections are read and we return as many cells as needed
% 
% References:
%   https://sites.google.com/site/silakovalexey/kazan-viewer/
%   A. Silakov Kazan viewer 2009
% See also: read_jeol, read_varian, read_bruker
  
  data = {};
  
  index = 1;
  while index > 0 && index < 6
      try
        this=opus_read(filename, index);
        if isempty(this), return;
        else data{end+1}=this; end
      end
      index = index +1;
  end
  
  

function ax=opus_read(varargin)
% OPUS FT-IR data file loader
%   [ax,y,ax]=opus_read(filename)
%   Loads one selected data section (see format description below )
% Made to be part of KAZAN Viewer
% alsi 24.11.2005
%
% FORMAT description
%   format is terrible. original data has no X-axis
%   Structure of the file is following:
%    first comes the difraction patern
%    second, FFT(baseline) cutted for some sertain region
%    and third FFT(data)-FFT(baseline)
%   If during experiment baseline is recordered (extention is usually *.0),
%   the third part is not present
%   Data is in the "float32" binary format. 
% HEADERS
%   Data sections are separated by headers, starting usually with "CSF" and end up
%   with "END". There is also some common header with description of
%   experimental setup Parameter, which I was able to find was "High/Low
%   folding limit" Headers have variable name in ascii code and values in
%   "float64" binary HFL/LFL - upper and lower frequency for X-axis, in
%   FFT(data) :-) 
%
% https://sites.google.com/site/silakovalexey/kazan-viewer/
% A. Silakov Kazan viewer 2009
%
%Copyright (c) 2009, A. Silakov
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%* Redistributions of source code must retain the above copyright
%  notice, this list of conditions and the following disclaimer.
%* Redistributions in binary form must reproduce the above copyright
%  notice, this list of conditions and the following disclaimer in
%  the documentation and/or other materials provided with the distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.


ax = [];
filename = varargin{1};
if nargin>1, datapart = varargin{2}; else datapart = 2; end
warning off MATLAB:nonIntegerTruncatedInConversionToChar;
fid=fopen(filename,'r', 'ieee-le');
if fid == -1
  return
end

%ntotal = 83236;
fseek(fid, 0, 'eof');
ntotal = ftell(fid);

frewind(fid);
if ntotal < 10e6
  tmpstr=char(fread(fid,ntotal,'schar')).'; % whole file in one string
else
  tmpstr=char(fread(fid,10e6,'schar')).'; % whole file in one string
end

% test if this is an OPUS file
hfl = strfind(tmpstr, 'HFL');
if isempty(hfl), fclose(fid); return; end
lfl = strfind(tmpstr, 'LFL');
if isempty(lfl), fclose(fid); return; end
fxv = strfind(tmpstr, 'FXV'); % first X point
lxv = strfind(tmpstr, 'LXV'); % last X point
npts = strfind(tmpstr, 'NPT'); % Number of data points 
if isempty(fxv) || isempty(lxv) || isempty(npts), fclose(fid); return; end

frewind(fid);
tmpdat=fread(fid,ntotal/4,'float32'); % whole file in one array

fseek(fid, hfl+7, 'bof');
HighFold = fread(fid, 1, 'float64');

fseek(fid, lfl+7, 'bof');
LowFold = fread(fid, 1, 'float64');


ndatas = min([length(fxv), length(lxv)]);
for ii = 1:ndatas
        fseek(fid, fxv(ii)+7, 'bof');    
    FreqStPoint(ii) = fread(fid, 1, 'float64');
        fseek(fid, lxv(ii)+7, 'bof');
    FreqEnPoint(ii) = fread(fid, 1, 'float64');
        fseek(fid, npts(ii)+7, 'bof');
    NDataPoints(ii) = fread(fid, 1, 'int32');
end

npts = strfind(tmpstr, 'RSN'); % Running sample number (size(data)?)
if ~isempty(npts)
    fseek(fid, npts+7, 'bof');
    NPoints = fread(fid, 1, 'int32');
else
    NPoints = 0;
end

%
%APT..2.5 mm 
%BMS KBr 
%CHN  Front 
%DTC  ...MCT 
%LPF  1  
%SRC MIR-Source 
%VEL  7   
%SGN -1  
%RGN -1  
%........END
%APF B3  
%HFQ  @?@
%LFQ �@ 
%PHR @
%PHZ ...ML  
%SPZ... NO  
%ZFF  2   
%........END     
%AQM  DN  
%COR  NO  
%DEL     
%DLY     
%HFW �@
%LFW     y@
%NSS   �  
%PLF AB  
%RES @ (bin 000001)
%TDL  &   
%.....END     
%CNM default 
% 
%EXP cryostat_mct4000-0_c.XPM ...Pun
%%% sample form
%SFM Cryostat 40K /2,5 mm  ..PM
%SNM Creinhardtii FeFeHase, HoxCO, rec, 200uM, 15uL ..K
%...END     
%CSF      �?
%MXY    ��?
%MNY `v�
%DAT 17/10/2007 
%TIM 03:09:36 ech
%DXU WN  
%.....END
npts = strfind(tmpstr, 'ARS'); % number of bg. scans
if ~isempty(npts) 
    fseek(fid, npts(1)+7, 'bof');
    NBGScans = fread(fid, 1, 'int32');
else
    NBGScans = 0;
end    
npts = strfind(tmpstr, 'ASS'); % number of sample scans
if ~isempty(npts) 
    fseek(fid, npts+7, 'bof'); 
    NScans = fread(fid, 1, 'int32');
else
    NScans = 0;
end
npts = strfind(tmpstr, 'PKL'); % Peak Location
if ~isempty(npts) 
    fseek(fid, npts+7, 'bof'); 
    PeakLocation = fread(fid, 1, 'int32');
else
    PeakLocation = 0;
end
npts = strfind(tmpstr, 'PRL'); % Backward Peak Location
if ~isempty(npts) 
    fseek(fid, npts+7, 'bof'); 
    BackwardPeakLocation = fread(fid, 1, 'int32');
else
    BackwardPeakLocation = 0;
end

% res = strfind(tmpstr, 'RES'); % Running sample number (size(data)?)
% fseek(fid, res+3, 'bof');
% Resol = fread(fid, 1, 'int32');

%aqm = strfind(tmpstr, 'AQM');

dd = strfind(tmpstr, 'DAT');
tt = strfind(tmpstr, 'TIM');
ax.Date = [tmpstr((dd(1)+8):(dd(1)+17))];
ax.Time = [tmpstr((tt(1)+8):(tt(1)+15))];

nptse = strfind(tmpstr, 'EXP'); 
nptsf = strfind(tmpstr, 'SFM'); 
nptsn = strfind(tmpstr, 'SNM'); 
nptsen  = strfind(tmpstr, 'END'); 
eie= find(nptsen>nptsn(1));
nptsen = nptsen(eie(1));
% nptsen = 
if ~isempty(nptsn)
    SampleName = [tmpstr((nptsn+8):(nptsen(1)-2))];
%     ee = isletter(SampleName);
%     iee = find(ee);
%     SampleName = SampleName(iee(1):iee(end));
else
    SampleName = '';
end
ax.SampleName = SampleName;
if ~isempty(nptse)
    ax.Experiment = [tmpstr((nptse+8):(nptsf-4))];
%     ee = isletter(ax.Experiment);
%     iee = find(ee);
%     ax.Experiment = ax.Experiment(iee(1):iee(end));
else
    ax.Experiment = '';
end
if ~isempty(nptsf)
    ax.SampleForm = [tmpstr((nptsf+8):(nptsn))];
%     ee = isletter(ax.SampleForm);
%     
%     iee = find(ee);
%     if ~isempty(iee)
%         ax.SampleForm = ax.SampleForm(iee(1):iee(end));
%     end
else
    ax.Experiment = '';
end

fclose(fid);

%stcom = strfind(tmpstr, 'CSF');
encom = strfind(tmpstr, 'END') + 3;

for ii = 1:length(fxv)
    tencom = encom(encom<fxv(ii));
    tencom = sort(tencom);
    if ~isempty(tencom)
      stcom(ii) = tencom(end);
    end
end
ax.Title = SampleName;
% 1-Interferogram, 2, 3 - FFT
if length(fxv)>2 % full set of datas
    datapart = min(length(fxv), datapart);
    if datapart==1,
        ax.xlabel = 'Mirror steps, pnts';
        ax.Title = [ax.Title, ' <S_IFG>'];
        ax.Label = 'SampleInterferogram';
    elseif datapart==2,
        ax.xlabel = 'Wavenumber, cm^{-1}';
        ax.Title = [ax.Title, ' <S_SC>'];
        ax.Label = 'SampleSpectrum';
    elseif datapart==3,
        ax.xlabel = 'Wavenumber, cm^{-1}';
        ax.Title = [ax.Title, ' <_AB>'];  
        ax.Label = 'RatioAbsorption';
    end
elseif length(fxv)==2 % baseline
    % somehow in file with baseline first comming fft(baseline) and only
    % then interferogram (sado_maso inc. :)
    datapart = min(length(fxv), datapart);
    if datapart==2,
        ax.xlabel = 'Wavenumber, cm^{-1}';
        ax.Title = [ax.Title, ' <S_SC>'];
        ax.Label = 'SampleSpectrum';
       datapart = 1;
    else
        ax.xlabel = 'Mirror steps, pnts';
        ax.Title = [ax.Title, ' <S_IFG>'];
        ax.Label = 'SampleInterferogram';
        datapart = 2;        
    end
elseif length(fxv)==1
    % only the AB
    datapart = min(length(fxv), datapart);
        ax.xlabel = 'Wavenumber, cm^{-1}';
        ax.Title = [ax.Title, ' <S_SC>'];
        ax.Label = 'SampleSpectrum';
end


y = tmpdat(floor((stcom(datapart)+3)/4)+1+(1 : NDataPoints(datapart)));

ax.x = linspace(FreqStPoint(datapart), FreqEnPoint(datapart), NDataPoints(datapart)).';
ax.y = y;
ax.Attributes.x = ax.xlabel;
ax.Attributes.y = ax.Label;

ax.HFL = HighFold;
ax.LFL = LowFold;

ax.HighFoldingLimit = [num2str(HighFold), ' cm-1'];
ax.LowFoldingLimit = [num2str(LowFold), ' cm-1'];
ax.NumberOfBGScans = [num2str(NBGScans)];
ax.NumberOfSampleScans = [num2str(NScans)];
ax.NumberOfDataSets = [num2str(length(fxv))];
ax.PeakLocation = num2str(PeakLocation);
ax.BackwardPeakLocation = num2str(BackwardPeakLocation);
