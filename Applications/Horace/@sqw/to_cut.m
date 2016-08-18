% Create cut object from sqw object
%
%   >> c = to_cut (w)
% With keywowrds:
%   >> c = to_cut (w, 'x','k')        % make x-axis the component along b* (other options available too)
%   >> c = to_cut (w, 'signal','Q')   % make y-axis equal to |Q| (other options available too)
%   >> c = to_cut (w, 'x', 'E', 'signal' 'k')
%
%   w           1D sqw object
%   'x'         [Optional] keyword to make the x-axis for the cut correspond to another coordinate
%               Valid coordinates are: 
%                   'h', 'k', 'l'       r.l.u.
%                   'E'                 energy transfer
%                   'Q'                 |Q|         
%               Default is to use the display axis of the sqw object
%   
%   'signal'    [Optional] keyword to make the signal axis another coordinate than the intensity.
%               Useful to see the variation of e.g. energy across a cut.
%   
%
% Note: this would normally be called just cut, which is the class of the output object.
% However, because we have already used 'cut' as the name of another method of sqw objects,
% we can't do that. An unfortunate problem, but one that is unavoidable.
%%   Overloaded methods:
%      sqw/to_cut
%      sqw/to_cut
%