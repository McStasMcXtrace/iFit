function s = mhdf4read(filename)
% mhdf4read Wrapper to hdfinfo/hdfread which reconstructs the HDF4 structure

s         = hdfinfo(filename);
if isfield(s, 'Raster8')
  s.Raster8 = hdfread(s.Raster8);
end
if isfield(s, 'Raster24')
  s.Raster24 = hdfread(s.Raster24);
end
if isfield(s, 'SDS')
  s.SDS = hdfread(s.SDS);
end
if isfield(s, 'Vdata')
  s.Vdata = hdfread(s.Vdata);
end
if isfield(s, 'Vgroup')
  s.Vgroup = hdfread(s.Vgroup);
end

