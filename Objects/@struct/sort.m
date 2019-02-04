function c = sort(s)
%  SORT Order fields of a structure array.
%  
%  SNEW = SORT(S1) orders the fields in S1 so the new structure array
%  SNEW has field names in ASCII dictionary order.

c = orderfields(s);
