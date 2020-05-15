function c = plus(a,b)
% +   Plus.
%   A + B adds matrices A and B.  
%   The addition is defined as the normalised sum:
%      (M1+M2)*(S1/M1+S2/M2) over monitor(M1+M2)
%   where S1,M1 and S2,M2 and the Signal and Monitor of the two objects.
%
%   To get the 'conventional' sum which is the sum S1+S2, set one of the
%   monitors to 0, e.g.:
%     a.Monitor=0;
%
%   Alternatively, the COMBINE operator (aka merge) of two objects is defined as:
%     (S1+S2) over monitor (M1+M2)
%
%   C = PLUS(A,B) is called for the syntax 'A + B'
%
% Example: a=estruct(peaks); b=a+2; max(b) == max(a)+2
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/minus, estruct/plus, estruct/times, estruct/rdivide, estruct/combine

if nargin ==1
	b=[];
end
c = binary(a, b, 'plus');

