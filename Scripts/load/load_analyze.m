function a = load_analyze(a)
%
% Returns an iData style dataset from an Analyze volume dataset
% typically used in medical imaging.
%
% Version: $Revision$
% See also: iData/load, iLoad, save, iData/saveas

if ~isa(a,'iData')
  a = load(iData,a,'Analyze');
  return
end
hdr=a.Data.hdr;
x = ([1:hdr.dim(1)]-hdr.origin(1))*hdr.siz(1);
y = ([1:hdr.dim(2)]-hdr.origin(2))*hdr.siz(2);
z = ([1:hdr.dim(3)]-hdr.origin(3))*hdr.siz(3);

setalias(a,'Signal','Data.img');
setalias(a,'x',x);
setalias(a,'y',y);
setalias(a,'z',z);
setaxis(a,1,'x');
setaxis(a,2,'y');
setaxis(a,3,'z');
