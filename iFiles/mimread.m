function s = mimread(filename)
% mimread Wrapper to imfinfo/imread which reconstructs the image structure

s      = imfinfo(filename);
s.data = imread(filename);

