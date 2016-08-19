function [u, v] = sw_cartesian(n)
% creates a right handed Cartesian coordinate system
%
% [u, v] = SW_CARTESIAN(n) 
%
% It creates an (n,u,v) right handed Cartesian coordinate system.
%
% n         3 element vector, either row or column.
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

if nargin == 0
    help sw_cartesian;
    return
end

% Shape of original vector.
nShape = size(n);

n = n(:);
z = [0; 0;-1];
y = [0;-1; 0];

if any(cross(n,z))
    u = cross(n,z);
else
    u = cross(n,y);
end

v = cross(n,u);
u = u/norm(u);
v = v/norm(v);

% Conserves the shape of the input vector.
u = reshape(u,nShape);
v = reshape(v,nShape);

end