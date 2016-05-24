function s = read_fits(filename)
% mfitsread Wrapper to fitsinfo/fitsread which reconstructs the FITS structure

try
  s      = fitsinfo(filename);
catch
  s=[]; return;
end
s.data = fitsread(filename);

