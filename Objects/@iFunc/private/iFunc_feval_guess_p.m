function [p,guessed,signal_in_varargin, ax] = iFunc_feval_guess_p(model, p, signal_in_varargin, varargin)

guessed = ''; ax = varargin;
model.Dimension = abs(model.Dimension);

% some ParameterValues have been defined already. Use them.
if isempty(p) && length(model.ParameterValues) == numel(model.Parameters)
  p = model.ParameterValues;
end
% when length(p) < Parameters, we fill NaN's ; when p=[] we guess them all
if strcmp(p, 'current'), p=model.ParameterValues; end
if isempty(p) % should guess parameters, but also evaluate model
  guessed = 'full and eval';
  p = NaN*ones(1, numel(model.Parameters));
elseif strcmp(p, 'guess') % explicitely return guessed parameters
  guessed = 'full';
  p = NaN*ones(1, numel(model.Parameters));
elseif isnumeric(p) && length(p) < length(model.Parameters) % fill NaN's from p+1 to model.Parameters
  if length(model.ParameterValues) == numel(model.Parameters)
    p((length(p)+1):length(model.Parameters)) = model.ParameterValues((length(p)+1):length(model.Parameters));
  else
    p((length(p)+1):length(model.Parameters)) = NaN;
  end
end

% when there are NaN values in parameter values, we replace them by guessed values
if model.Dimension && ...
  ((any(isnan(p)) && length(p) == length(model.Parameters)) || ~isempty(guessed))
  % call private method to guess parameters from axes, signal and parameter names
  if isempty(guessed), guessed = 'partial'; end
  
  % args={x,y,z, ... signal}
  args=cell(1,model.Dimension+1); args(1:end) = { [] };
  args(1:min(length(ax),model.Dimension+1)) = ax(1:min(length(ax),model.Dimension+1));
  args_opt = ax((model.Dimension+2):end);
  
  p0 = p; % save initial 'p' values
  
  % all args are empty, we generate model fake 1D/2D axes/signal
  if all(cellfun('isempty',args))
    if model.Dimension == 1 % we use a Gaussian
      args{1} = linspace(-5,5,50); 
      x=args{1}; p2 = [1 mean(x) std(x)/2 .1]; 
      args{2} = p2(1)*exp(-0.5*((x-p2(2))/p2(3)).^2)+((p2(2)-x)*p2(1)/p2(3)/100) + p2(4);
      clear p2
      signal = args{2};
      signal_in_varargin = 2;
    elseif model.Dimension == 2 % we use a 2D Gaussian
      [args{1},args{2}] = ndgrid(linspace(-5,5,30), linspace(-3,7,40));
      x=args{1}; y=args{2}; p2 = [ 1 mean(x(:)) mean(y(:)) std(x(:)) std(y(:)) 30 0 ];
      x0=p2(2); y0=p2(3); sx=p2(4); sy=p2(5);
      theta = p2(6)*pi/180;  % deg -> rad
      aa = cos(theta)^2/2/sx/sx + sin(theta)^2/2/sy/sy;
      bb =-sin(2*theta)/4/sx/sx + sin(2*theta)/4/sy/sy;
      cc = sin(theta)^2/2/sx/sx + cos(theta)^2/2/sy/sy;
      args{3} = p2(1)*exp(-(aa*(x-x0).^2+2*bb*(x-x0).*(y-y0)+cc*(y-y0).^2)) + p2(7);
      clear aa bb cc theta x0 y0 sx sy p2
      signal = args{3};
      signal_in_varargin = 3;
    else % use an event style representation
      for index=1:(model.Dimension+1)
        x1 = -2*rand-1;
        x2 = 2*rand+1;
        args{index} = linspace(x1, x2, 20+index);
      end
      signal_in_varargin = model.Dimension+1;
      signal = args{end};
    end
  end
  
  ax = [ args args_opt ];
  clear args
  
  % convert axes to nD arrays for operations to take place
  % check the axes and possibly use ndgrid to allow nD operations in the
  % Expression/Constraint
  % Not for event style axes+signal (all 1D)
  
  % event:  all vectors, including signal (if any), same length
  % regrid: all vectors, not same length, signal is not a vector
  
  % check for vectors: 0 or 1 for each axis
  check_vector = cellfun(@(c)(~isempty(c) && length(c) == numel(c)), ax(1:model.Dimension));
  % check number of elements: would be equal for 3D grids (xyz) or 4D (xyzt)
  check_numel  = cellfun(@numel, ax(1:model.Dimension));
  % check orientation (axis rank) for vectors. 0 for non vectors
  check_orient = zeros(1,model.Dimension);
  index        = find(check_vector);
  try
    check_orient(index) = cellfun(@(c)find(size(c)==numel(c)), ax(index));
  catch
    check_orient(index) = 0;
  end
  
  % all vectors, not same orientation and not same length: ndgrid
  if all(check_vector) && any(check_orient ~= check_orient(1)) ...
    && ~all(check_numel  == check_numel(1) | check_numel == 1)
    [ax{1:model.Dimension}] = ndgrid(ax{1:model.Dimension});
  end
  % automatic guessed parameter values -> signal
  if model.Dimension
    p1 = iFunc_private_guess(ax(1:(model.Dimension+1)), model.Parameters); % call private here -> auto guess
  else
    p1 = [];
  end
  % check for NaN guessed values and null amplitude
  n=find(isnan(p1) | p1 == 0); n=transpose(n(:));
  for j=n
    if any([strfind(lower(model.Parameters{j}), 'width') ...
       strfind(lower(model.Parameters{j}), 'amplitude') ...
       strfind(lower(model.Parameters{j}), 'intensity')])
      p1(j) = abs(randn)/10;
    else
      p1(j) = 0;
    end
  end
  % specific guessed values (if any) -> p2 override p1
  if ~isempty(model.Guess) && ~all(cellfun('isempty',ax))
    if ischar(model.Guess)
      % request char eval guess in sandbox
      p2 = iFunc_feval_guess(model, ax{:});
      if isa(p2, 'function_handle')
        model.Guess = p2;
      end
    end
    if isa(model.Guess, 'function_handle')
      try
        n = nargin(model.Guess);                % number of required arguments
        % moments of distributions
        m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
        m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));
        if n == ndims(model)+2
            % syntax Guess: @(p,x,y,... signal)
            if n > 0 && length(ax) >= n
              p2 = feval(model.Guess, p, ax{1:n+1}); % returns model vector
            else
              p2 = feval(model.Guess, p, ax{:}); % returns model vector
            end
        else
            % syntax Guess: @(x,y,... signal)
          if n > 0 && length(ax) >= n
            p2 = feval(model.Guess, ax{1:n}); % returns model vector
          else
            p2 = feval(model.Guess, ax{:}); % returns model vector
          end
        end
      catch ME
        disp([ mfilename ': Guess: ' ME.message ])
        p2 = [];
      end
      clear n
    elseif isnumeric(model.Guess)
      p2 = model.Guess;
    else
      p  = p0;             % restore initial value
    end
    if isempty(p2)
      disp([ mfilename ': Warning: Could not evaluate Guess in model ' model.Name ' ' model.Tag ]);
      disp(model.Guess);
      disp('Axes and signal:');
      disp(ax);
      warning('Using auto-guess values.');
    else
      % merge auto and possibly manually set values
      index     = ~isnan(p2);
      p1(index) = p2(index);
      clear p2
    end
  end
  if all(p1 == 0) && ~isempty(model.ParameterValues) ...
   && ~all(model.ParameterValues(:) == 0)
    p1 = model.ParameterValues;
  end
  signal = p1;  % auto-guess overridden by 'Guess' definition
  % transfer the guessed values from 'signal' to the NaN ones in 'p'
  if any(isnan(p)) && ~isempty(signal)
    index = find(isnan(p)); p(index) = signal(index);
  end
  model.ParameterValues = p; % the guessed values
  
  if ~strcmp(guessed,'full')
    % return the signal and axes
    % [signal, ax, name] = feval(model, p, varargin{:});
    guessed = ''; % proceed with eval
  end
  % Parameters are stored in the updated model (see assignin below)
end % 'guess'
