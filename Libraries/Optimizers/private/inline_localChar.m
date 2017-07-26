function strfcn = inline_localChar(fcn)
% Convert the fcn to a string for printing

  if ischar(fcn)
      strfcn = fcn;
  elseif isa(fcn,'inline')
      strfcn = char(fcn);
  elseif isa(fcn,'function_handle')
      strfcn = func2str(fcn);
  else
      try
          strfcn = char(fcn);
      catch
          strfcn = '(name not printable)';
      end
  end

end % inline_localChar
