function out = iData_check(in)
% iData_chec: make consistency checks on iData object

if numel(in) > 1
  out = zeros(iData, numel(in), 1);
  parfor index = 1:numel(in)
    out(index) = iData_check(in(index));
  end
  return
end

if iscell(in), in = in{1}; end

% update ModifDate
in.ModificationDate = clock;
% check type of fields
if iscellstr(in.Title)
  t = strcat(in.Title,';');
  in.Title=[ t{:} ];
end
if ~ischar(in.Title) 
  iData_private_warning(mfilename,['Title must be a char or cellstr in iData object ' in.Tag ' (' class(in.Title) '). Re-setting to empty.']);
  in.Title = '';
end
in.Title = strtrim(in.Title);
if ~ischar(in.Tag)
  iData_private_warning(mfilename,['Tag must be a char in iData object ' in.Tag ' "' in.Title '. Re-setting to a new Tad id.' ]);
  in = iData_private_newtag(in);
end
if ~ischar(in.Source)
  iData_private_warning(mfilename,['Source must be a char in iData object ' in.Tag ' "' in.Title '. Re-setting to pwd.' ]);
  in.Source = pwd;
end
if ~iscellstr(in.Command)
  in.Command = { in.Command };
end
if ~ischar(in.Date) && ~isnumeric(in.Date)
  iData_private_warning(mfilename,['Date must be a char/vector in iData object ' in.Tag ' "' in.Title '. Re-setting to "now".' ]);
  in.Date = clock;
end
if ~ischar(in.Creator)
  iData_private_warning(mfilename,['Creator must be a char in iData object ' in.Tag ' "' in.Title '. Re-setting to "iFit/iData".' ]);
  in.Creator = version(in);
end
if ~ischar(in.User)
  iData_private_warning(mfilename,['User must be a char in iData object ' in.Tag ' "' in.Title '. Re-setting to Matlab User.']);
  in.User = 'Matlab User';
end
% check if object.Data is numeric: make it a structure so that it is better organized
if isnumeric(in.Data) && ~isempty(in.Data)
  data = in.Data;
  in.Data = [];
  in.Data.Signal = data;
end

if ~isempty(in.Data) && isempty(getalias(in, 'Signal'))
  % get numeric fields sorted in descending size order
  [fields, types, dims] = findfield(in, '', 'numeric');
  if isempty(fields), 
    iData_private_warning(mfilename,['The iData object ' in.Tag ' "' in.Title '" contains no data at all ! (double/single/logical/int/uint)']);
  else
    fields_all = fields; dims_all=dims;
    % does this looks like a Signal ?
    
    if length(dims) > 1 % when similar sizes are encoutered, get the one which is not monotonic
      % list of 'biggest' fields
      maxdim=find(dims == dims(1)); maxdim2 = maxdim;
      % move 'error' and constant/monotonic down in the list
      for idx=1:length(maxdim)
        index=maxdim2(idx);
        x = get(in, fields{index});
        if ischar(x) || length(x) <= 1
          % this is a char: move to end of list
          maxdim([ end idx]) = maxdim([ idx end] );
          continue; 
        end
        x = diff(x(:));
        if all(x == x(1)) || all(x > 0) ...
          || ~isempty(strfind(lower(fields{index}), 'error')) ...
          || ~isempty(strfind(lower(fields{index}), 'monitor'))
          % this is a constant/monotonic value or 'error' or 'monitor'
          % move down in fields list
          maxdim([ end idx]) = maxdim([ idx end] );
        end
      end
      fields(maxdim2) = fields(maxdim);
    end
    
    % in case we have more than one choice, get the first one and error bars
    error_id = []; monitor_id=[];
    if length(dims) > 1 || iscell(fields)
      % do we have an 'error' which has same dimension ?
      for index=find(dims(:)' == dims(1))
        if index==1, continue; end % not the signal itself
        if ~isempty(strfind(lower(fields{index}), 'error'))
          error_id = fields{index};
        elseif ~isempty(strfind(lower(fields{index}), 'monitor'))
          monitor_id = fields{index};
        end
      end
      
      dims=dims(1);
      fields=fields{1};
    end

    % index: is the field and dimension index to assign the Signal
    if dims > 0
      disp([ 'iData: Setting Signal="' fields '" with length ' num2str(dims) ' in object ' in.Tag ' "' in.Title '".' ]);
      in = setalias(in,'Signal', fields);
      
      % get potential attribute (in Data.Headers or Data.Attributes fields)
      attribute = iData_getAttribute(in, fields);
      
      if isstruct(attribute)
        attribute = class2str(' ',attribute, 'no comments');
      end
      if ischar(attribute)
        if numel(attribute) > 80, attribute=[ attribute(1:77) ' ...' ]; end
        in.Alias.Labels{1} = attribute;
      end

      % assign potential 'error' bars and 'monitor'
      if ~isempty(error_id)
        in = setalias(in,'Error', error_id);
      end
      if ~isempty(monitor_id)
        in = setalias(in,'Monitor', monitor_id);
      end
    end
    % look for vectors that may have the proper length as axes
    for index=1:ndims(in)
      if isempty(getaxis(in, num2str(index)))
        % search for a vector of length size(in, index)
        ax = find(dims_all == size(in, index));   % length of dim, or length(dim)+1
        if isempty(ax), ax = find(dims_all == size(in, index)+1); end
        ax = ax(~strcmp(getalias(in,'Signal'), fields_all(ax)));
        if length(ax) > 1; ax=ax(1); end
        if ~isempty(ax)
          val = get(in, fields_all{ax});
          if isvector(val) && ~strcmp(fields_all{ax},getalias(in,'Signal'))
            if length(val) == size(in, index) && min(val(:)) < max(val(:))
              in = setaxis(in, index, [ 'Axis_' num2str(index) ], fields_all{ax});
              found = 1;
            elseif length(val) == size(in, index)+1 && min(val(:)) < max(val(:))
              val = (val(1:(end-1)) + val(2:end))/2;
              in = setaxis(in, index, [ 'Axis_' num2str(index) ], val);
              found = 1;
            else found = 0;
            end
            if found == 1  % the axis could be found
              % search if there is a corresponding label (in Headers)
              if isfield(in.Data, 'Headers')
                fields=fliplr(strtok(fliplr(fields_all{ax}), '.'));
                if isfield(in.Data.Headers, fields)
                  in.Alias.Labels{index+1} = in.Data.Headers.(fields);
                end
              else
                label(in, index, fields_all{ax});
              end
              disp([ 'iData: Setting Axis{' num2str(index) '} ="' fields_all{ax} '" with length ' num2str(length(val)) ' in object ' in.Tag ' "' in.Title '".' ]);
            end
          end
          clear val
        else
          break; % all previous axes must be defined. If one misses, we end the search
        end
      end
    end 
  end
end

% check in case the x,y axes have been reversed for dim>=2, then swap 1:2 axes
if ndims(in)==2 && ~isempty(getaxis(in, '1')) && ~isempty(getaxis(in, '2')) ...
            && isvector(getaxis(in, 1)) && isvector(getaxis(in, 2)) ...
            && length(getaxis(in, 1)) == size(get(in,'Signal'),2) ...
            && length(getaxis(in, 2)) == size(get(in,'Signal'),1) ...
            && length(getaxis(in, 1)) ~= length(getaxis(in, 2))
  x1 = getaxis(in, '1');
  x2 = getaxis(in, '2');
  setaxis(in, 1, x2);
  setaxis(in, 2, x1);
  clear x1 x2
  disp([ 'iData: The object has been transposed to match the axes orientation in object ' in.Tag ' "' in.Title '".' ]);
end

% check aliases (valid ?) by calling setalias(in)
in = setalias(in);
% check axis (valid ?) by calling setaxis(in)
in = setaxis(in);

out = in;


