% Create a 0D Horace dataset ('d0d')
%
% Syntax:
%   >> w = d0d (filename)       % Create object from a file
%
%   >> w = d0d (din)            % Create from a structure with valid fields
%                               % Structure array will output an array of objects
%
% Or:
%   >> w = d0d (u0)             % u0 is offset of origin of dataset,
%   >> w = d0d (lattice, u0)    % Give lattice parameters (a,b,c,alf,bet,gam)
%
% Input parameters in more detail:
% ----------------------------------
%   lattice Defines crystal lattice: [a,b,c,alpha,beta,gamma]
%   u0      Vector of form [h0,k0,l0] or [h0,k0,l0,en0]
%          that defines an origin point on the manifold of the dataset.
%          If en0 omitted, then assumed to be zero.
%