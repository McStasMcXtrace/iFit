function this = private_postprocess(this, postprocess)
% evaluate the postprocess in a reduced environment, with only 'this'
  this0 = this;
  if isempty(postprocess), return; end
  try
    if isvarname(postprocess) && exist(postprocess) == 2
      this = feval(postprocess, this);
    elseif ~any(postprocess == '=')
      this = eval(postprocess);
    else
      eval(postprocess);
    end
    
    % restore initial object 'this' if it was destroyed at eval
    if ~isa(this, 'estruct'), this = this0; end
    this = history(this, postprocess, this);

  catch ME
    if this0.verbose
      warning(getReport(ME));
      warning([mfilename ': Error when calling post-process ' postprocess '. file: ' this.Source ]);
    end
  end
