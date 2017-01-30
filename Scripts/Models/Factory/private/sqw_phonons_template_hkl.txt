function [HKL, t, is_event, resize_me] = sqw_phonons_template_hkl(x,y,z,t)
% sqw_phonons_template_hkl: code to analyze [xyzt] and prepare a list of HKL locations
% can work in the following logic:
% event 4D: [ x y z t ]   as vectors, same length, same orientation -> ResLibCal...
% kpath 2D: [ x y z ]     as vectors, same length, same orientation, [ t ] not same length or orientation
% map 3D:   [ x y z ]     are all 3D arrays, and [ t ] gives 4D as a vector
% map 4D:   any other choice. Get unique values, and build ndgrid

ax = {x y z t};
% check for vectors: 0 or 1 for each axis
check_vector = cellfun(@(c)(~isempty(c) && length(c) == numel(c)), ax);
% check number of elements: would be equal for 3D grids (xyz) or 4D (xyzt)
check_numel  = cellfun(@numel, ax);
% check orientation (axis rank) for vectors. 0 for non vectors
check_orient = zeros(1,4);
index        = find(check_vector);
check_orient(index) = cellfun(@(c)find(size(c)==numel(c)), ax(index));

HKL = []; is_event = false; resize_me = [];
if all(check_vector) && all(check_numel == check_numel(1)) && all(check_orient == check_orient(1))
  % event 4D: [ x y z t ]   as vectors, same length, same orientation
  is_event = true;
  
elseif all(check_vector(1:4)) && all(check_numel(1:3) == check_numel(1)) && all(check_orient(1:3) == check_orient(1))
  % kpath 2D: [ x y z ]     as vectors, same length, same orientation, [ t ] vector, not same length or orientation
  % this is the SAME test as event 4D, but with xyz only

elseif all(cellfun(@ndims, ax) == 4) && all(check_numel == check_numel(1)) % (X,Y,Z,T) are all 4D, same numel
  x=x(:,:,:,1); y=y(:,:,:,1); z=z(:,:,:,1); t=t(1,1,1,:);
  resize_me = size(x);

elseif all(cellfun(@ndims, ax(1:3)) == 3)  && all(check_numel(1:3) == check_numel(1)) % (X,Y,Z) are all 3D
  % 3D HKL data sets plus energy axis. NOP.
  resize_me = size(x);
else
  % map 4D:   any other choice
  % any(~check_vector) || any(check_numel ~= check_numel(1)) || any(check_orient ~= check_orient(1))
  x=unique(x); y=unique(y); z=unique(z); 
  [x,y,z] = ndgrid(x,y,z);
  resize_me = size(x);
end
clear ax
HKL      = [ x(:) y(:) z(:) ];      % a clean array (Nx3) of HKL coordinates
diff_HKL = diff(HKL,1,1);           % difference along rows
index = find(sum(diff_HKL.^2,2) == 0);
HKL(index, :) = [];  % remove consecutive duplicates at same HKL
x(index)= []; y(index)= []; z(index) = [];
t        = unique(t);
clear diff_HKL
if numel(resize_me) == 3
  resize_me = [ resize_me numel(t) ];
end

% now we compute FREQ (and POLAR when possible) at HKL
POLAR=[]; T=0; Amplitude=1; Gamma=0.01;  % defaults
