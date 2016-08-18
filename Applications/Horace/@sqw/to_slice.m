% Create slice object from sqw object
%
%   >> s = to_slice (w)
%   >> s = to_slice (w, 'signal','Q')   % make signal equal to |Q| (other options available too)
%
%   w           2D sqw object
%   
%   'signal'    [Optional] keyword to make the signal axis another coordinate than the intensity.
%               Useful to see the variation of e.g. energy across a slice.
%               Valid coordinates are: 
%                   'h', 'k', 'l'       r.l.u.
%                   'E'                 energy transfer
%                   'Q'                 |Q|         
%%   Overloaded methods:
%      sqw/to_slice
%      sqw/to_slice
%