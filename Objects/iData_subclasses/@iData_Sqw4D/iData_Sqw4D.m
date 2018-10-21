classdef iData_Sqw4D < iData
  % iData_Sqw4D: create a 4D S(q,w) data set (iData flavour)
  %
  % The iData_Sqw4D class is a 4D data set holding a S(h,k,l,w) dynamic structure factor
  %   aka scattering function/law.
  %   The axes are
  %   1 QH (rows)    Momentum transfer along H [rlu] (wavevector). 
  %   2 QK (columns) Momentum transfer along K [rlu] (wavevector). 
  %   3 QK (pages)   Momentum transfer along K [rlu] (wavevector). 
  %   4 EN           Energy transfer [meV]. 
  %
  % conventions:
  % w = omega = Ei-Ef = energy lost by the neutron [meV]
  %    omega > 0, neutron looses energy, can not be higher than Ei (Stokes)
  %    omega < 0, neutron gains energy, anti-Stokes
  %
  % Useful methods for this iData flavour:
  %
  % methods(iData_Sqw4D)
  %   all iData methods can be used.
  % iData_Sqw4D(s)
  %   convert input [e.g. a 4D iData object] into an iData_Sqw4D to give access to
  %   the methods below.
  % p   = parseparams(s)
  %   Search for physical quantities in S(q,w) data set.
  % saveas(s, filename, 'McStas')
  %   Save the S(q,w) as a McStas Sqw 4D, or other file format
  %
  % See also: iData, iData_Sqw2D, iData_vDOS, iFunc_Sqw4D
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
      if ischar(s), s = iData(s); end
      % convert/test
      if     isa(s, mfilename)   m = s;
      else
        if ~isa(s, 'iData') || any(isempty(s)) || any(ndims(s) ~= 4)
          error([ mfilename ': the given input ' class(s) ' does not seem to be convertible to iData_Sqw4D.' ])
        end
      end
      
      % copy all properties
      obj = copy_prop(obj, s);
      obj = commandhistory(obj, mfilename, s);
      label(obj, 0, [  mfilename '(' label(obj, 0) ')' ]);
      
      % search for Sqw parameters for further conversions
      [obj, parameters] = Sqw_parameters(obj);
 
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
    
    % parameters (search for parameters in iData)
    function parameters = parseparams(s)
      % iData_Sqw4D: parseparams: search for physical quantities in object.
      % This search is also done when creating iData_Sqw4D objects.
      %
      %   iData_Sqw4D -> physical parameters
      [s,parameters,fields] = Sqw_parameters(s);
      if length(inputname(1))
        assignin('caller',inputname(1),s);
      end
    end
    
    function f = saveas(self, varargin)
      % iData_Sqw4D: saveas: save S(h,k,l,w) into a file.
      %
      % convert: iData_Sqw4D S(h,k,l,w) -> file
      %
      % syntax: saveas(sqw4D, filename, format)
      %
      %   with format as mcstas, sqw or any other iData supported format.
      %
      % all iData.saveas formats are available, and in addition:
      %
      %   McStas (for Single_crystal_inelastic component)
      %
      if numel(varargin) >= 2
        switch lower(varargin{2})
        case {'mcstas','sqw'}
          f = Sqw_McStas(self, varargin{1});
        otherwise
          f = saveas@iData(self, varargin{:});
        end
      else
        f = saveas@iData(self, varargin{:});
      end
      
    end
    
  end
  
end
