function w = d0d (varargin)
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


% Original author: T.G.Perring
%
% $Revision: 301 $ ($Date: 2009-11-03 20:52:59 +0000 (Tue, 03 Nov 2009) $)


ndim_request = 0;
class_type = 'd0d';
inferiorto('sqw');

% ----- The following shoudld be independent of d0d, d1d,...d4d ------------
% Work via sqw class type

if nargin==1 && isa(varargin{1},class_type)     % already object of class
    w = varargin{1};
    return
end

if nargin==0
    w = sqw('$dnd',ndim_request); % default constructor
else
    w = sqw('$dnd',varargin{:});
    if dimensions(w)~=ndim_request
        error(['Input arguments inconsistent with requested dimensionality ',num2str(ndim_request)])
    end
end

w=class(w.data,class_type);
