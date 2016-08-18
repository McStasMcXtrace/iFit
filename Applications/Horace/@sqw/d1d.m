% Create a 1D Horace dataset ('d1d')
%
% Syntax:
%   >> w = d1d (filename)       % Create object from a file
%
%   >> w = d1d (din)            % Create from a structure with valid fields
%                               % Structure array will output an array of objects
%
% Or:
%   >> w = d1d (u1,p1      )    % u1,u2 vectors define projection axes in rlu,
%                                 p1,p2 give start,step and finish for the axes
%   >> w = d1d (u0,...)         % u0 is offset of origin of dataset,
%   >> w = d1d (lattice,...)    % Give lattice parameters [a,b,c,alf,bet,gam]
%   >> w = d1d (lattice,u0,...) % Give u0 and lattice parameters
%
% Input parameters in more detail:
% ----------------------------------
%   lattice Defines crystal lattice: [a,b,c,alpha,beta,gamma]
%   u0      Vector of form [h0,k0,l0] or [h0,k0,l0,en0]
%          that defines an origin point on the manifold of the dataset.
%          If en0 omitted, then assumed to be zero.
%   u1      Vector [h1,k1,l1] or [h1,k1,l1,en1] defining a plot axis. Must
%          not mix momentum and energy components e.g. [1,1,2], [0,2,0,0] and
%          [0,0,0,1] are valid; [1,0,0,1] is not.
%   p1      Vector of form [plo,delta_p,phi] that defines limits and step
%          in multiples of u1.
%