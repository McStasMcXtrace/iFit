classdef iData_Sqw4D < iData
  % iData_Sqw4D: create a 4D S(q,w) data set (iData flavour)
  %
  % The iData_Sqw4D class is a 4D data set holding a S(h,k,l,w) dynamic structure factor
  %   aka scattering function/law.
  %
  % conventions:
  % w = omega = Ei-Ef = energy lost by the neutron [meV]
  %    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
  %    omega < 0, neutron gains energy, anti-Stokes
  %
  % See also: iData, iData_Sab, iData_vDOS
  % (c) E.Farhi, ILL. License: EUPL.

  properties
  end
  
  methods
    % main instantiation method
    function obj = iData_Sqw4D(s)
      % iData_Sqw4D: create the iData_Sqw4D subclass
      %
      %   convert: anything -> iData_Sqw4D
      %
      % input:
      %   must be a 4D iData object holding a S(h,k,l,w)
      %   or file name.
      %
      % Example: s=iData_Sqw4D(file);
      obj = obj@iData;
      obj.class = mfilename;
      
      if ~nargin, return; end  % empty object
      
      % convert/test
      if     isa(s, mfilename)   m = s;
      else
        if ~isa(m, 'iData') || any(isempty(m)) || any(ndims(m) ~= 4)
          error([ mfilename ': the given input ' class(s) ' does not seem to be convertible to iData_Sqw4D.' ])
        end
      end
      
      % copy all properties
      obj = copy_prop(obj, m);
      obj = commandhistory(obj, mfilename, s);
      label(obj, 0, [  mfilename '(' label(obj, 0) ')' ]);
 
    end % iData_Sqw4D constructor
    
    function f = iData(self)
      % iData_Sqw4D: iData: convert a iData_Sqw4D back to iData
      %
      % convert: iData_Sqw4D -> iData
      f = [];
      for index=1:numel(self)
        f1   = copy_prop(iData, self(index));
        f1   = commandhistory(f1, 'iData', self(index));
        label(f1, 0, [  'iData' '(' label(self(index), 0) ')' ]);
        f = [ f f1 ];
      end
    end
    
  end
  
end
