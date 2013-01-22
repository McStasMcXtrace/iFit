function s = read_image(filename)
% mimread Wrapper to imfinfo/imread which reconstructs the image structure

s       = imfinfo(filename);
s.image = imread(filename);
try
  s.EXIF = exifread(filename);
end

