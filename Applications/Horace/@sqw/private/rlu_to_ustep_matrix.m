%  Perform various calculations with reciprocal lattice for producing plot axes
%
% Syntax:
%   >> [rlu_to_ustep, u_to_rlu, ulen] = rlu_to_ustep_matrix (alatt, angdeg, u, v, ustep, type, w)
%
% input:
% -----------
%   alatt(1:3)  Row vector of lattice parameters (Angstroms)
%   angdeg(1:3) Row vector of lattice angles (degrees)
%   u(1:3)      Row vector defining first plot axis (r.l.u.)
%   v(1:3)      Row vector defining plane of plot in Q-space (r.l.u.)
%           The plot plane is defined by u and the perpendicular to u in the
%           plane of u and v. The unit lengths of the axes are determined by the
%           character codes in the variable 'type' described below
%            - if 'a': unit length is one inverse Angstrom
%            - if 'r': then if (h,k,l) in r.l.u., is normalised so max(abs(h,k,l))=1
%            - if 'p': then normalised so that if the orthogonal set created from u and v is u1, u2, u3:
%                       |u1|=|u|, (u x u2)=(u x v), (u x u3)=(u x w)
%                      i.e. the projections of u,v,w along u1,u2,u3 match the lengths of u1,u2,u3
%
%   ustep(1:3)  Row vector giving step size along u1, u2 and u3 axes
%   type        Units of binning and thickness: a three-character string,
%              each character indicating if u1, u2, u3 normalised to Angstrom^-1
%              or r.l.u., max(abs(h,k,l))=1 - 'a' and 'r' respectively. e.g. type='arr'
%   w(1:3)      Row vector defining the line of the third axis. Only needed if type(3)='p' (r.l.u.)
%
% output:
% -----------
%   rlu_to_ustep(1:3,1:3)   Matrix to convert components of a vector expressed
%                           in r.l.u. to multiples of the step size along the
%                           orthogonal set defined by the vectors u and v
%                       i.e.
%                           Vstep(i) = rlu_to_ustep(i,j)*Vrlu(j)
%
%   u_to_rlu(1:3,1:3)       Vectors u1, u2, u3 in reciprocal lattice vectors: 
%                           the ith column is ui i.e. ui = u_to_rlu(:,i) 
%
%   ulen(1:3)               Row vector of lengths of ui in Ang^-1
%
%   mess                    Error message
%                           - all OK:   empty
%                           - if error: message, and rlu_to_ustep, u_to_rlu, ulen are empty
%                           
%