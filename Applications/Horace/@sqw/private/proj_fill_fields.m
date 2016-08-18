% Fill the fields of proj data sstructure from input and defaults for missing information
%
% Input:
%   proj_in Data structure containing details of projection axes:
%           Defines two vectors u and v that give the direction of u1
%           (parallel to u) and u2 (in the plane of u1 and u2, with u2
%           having positive component along v); also defines the 
%           normalisation of u1,u2,u3
%             Required arguments:
%               proj.u          [1x3] Vector of first axis (r.l.u.)
%               proj.v          [1x3] Vector of second axis (r.l.u.)
%             Optional arguments:
%               proj.w          [1x3] Vector of third axis (r.l.u.) - only needed if third character of type is 'p'
%                               Will otherwise be ignored.
%               proj.type       [1x3] Char. string defining normalisation:
%                   Each character indicates if u1, u2, u3 normalised as follows:
%                   - if 'a': unit length is one inverse Angstrom
%                   - if 'r': then if (h,k,l) in r.l.u., is normalised so max(abs(h,k,l))=1
%                   - if 'p': then normalised so that if the orthogonal set created from u and v is u1, u2, u3:
%                       |u1|=|u|, (u x u2)=(u x v), (u x u3)=(u x w)
%                      i.e. the projections of u,v,w along u1,u2,u3 match the lengths of u1,u2,u3
%                   Default:
%                       'ppr'  if w not given
%                       'ppp'  if w is given
%               proj.uoffset    Row or column vector of offset of origin of projection axes (r.l.u.)
%               proj.lab        Short labels for u1,u2,u3,u4 as cell array (e.g. {'Q_h', 'Q_k', 'Q_l', 'En'})
%                  *OR*
%               proj.lab1       Short label for u1 axis (e.g. 'Q_h' or 'Q_{kk}')
%               proj.lab2       Short label for u2 axis
%               proj.lab3       Short label for u3 axis
%               proj.lab4       Short label for u4 axis (e.g. 'E' or 'En')
%
% Output:
%   proj    Output data structure with defaults for absent fields (note that labels packed as single cellstr)
%              proj.u          [1x3] Vector of first axis (r.l.u.)
%              proj.v          [1x3] Vector of second axis (r.l.u.)
%           Optional arguments:
%              proj.w          [1x3] Vector of third axis (r.l.u.) (set to [] if not given in proj_in)
%              proj.type       [1x3] Char. string defining normalisation:
%              proj.uoffset    [4x1] column vector of offset of origin of projection axes (r.l.u. and en)
%              proj.lab        [1x4] cell array of projection axis labels
%