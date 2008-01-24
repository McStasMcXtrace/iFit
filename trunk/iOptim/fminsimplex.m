function [pars,fval,exitflag,output] = fminsimplex(fun, pars, options)
% [MINIMUM,FVAL,EXITFLAG,OUTPUT] = FMINSIMPLEX(FUN,PARS,[OPTIONS]) Nelder-Mead Simplex
%
% This minimization method uses an adaptive random search.
% 
% Calling:
%   fminsimplex(fun, pars) asks to minimize the 'fun' objective function with starting
%     parameters 'pars' (vector)
%   fminsimplex(fun, pars, options) same as above, with customized options (optimset)
%
% Example:
%   banana = @(x)100*(x(2)-x(1)^2)^2+(1-x(1))^2;
%   [x,fval] = fminsimplex(banana,[-1.2, 1])
%
% Input:
%  FUN is the function to minimize (handle or string).
%
%  PARS is a vector with initial guess parameters. You must input an
%  initial guess.
%
%  OPTIONS is a structure with settings for the optimizer, 
%  compliant with optimset. Default options may be obtained with
%   optimset('fminsimplex')
%
% Output:
%          MINIMUM is the solution which generated the smallest encountered
%            value when input into FUN.
%          FVAL is the value of the FUN function evaluated at MINIMUM.
%          EXITFLAG return state of the optimizer
%          OUTPUT additional information returned as a structure.
% Reference: Nelder and Mead, Computer J., 7 (1965) 308
% Contrib: F. Sigworth, 15 March 2003, S. H. Heinemann, 1987
%          M. Caceci and W. Cacheris, Byte, p. 340, May 1984.
%
% Version: $Revision: 1.3 $
% See also: fminsearch, optimset

% default options for optimset
if nargin == 1 & strcmp(fun,'defaults')
  options=optimset; % empty structure
  options.Display='off';
  options.TolFun =1e-4;
  options.TolX   =1e-6;
  options.MaxIter=500;
  options.MaxFunEvals=501;
  pars = options;
  return
end

if nargin <= 2
  options=[];
end
if isempty(options)
  options=feval(mfilename, 'defaults');
end

if strcmp(options.Display,'iter')
  fmin_private_disp_start(mfilename, fun, pars);
end

[p,t]=Simplex('init',pars);
iters=0; funcount=0;
while 1
  p_prev=p;
  y=feval(fun, p);
  funcount=funcount+1;
  iters=iters+1;
  [p,t]=Simplex(y);
  % std stopping conditions
  [istop, message] = fmin_private_std_check(p, y, iters, funcount, options, p_prev);
  if strcmp(options.Display, 'iter')
    fmin_private_disp_iter(iters, funcount, fun, p, y);
  end
  if Simplex('converged',options.TolFun)
    istop=-9;
    message = [ 'Global Simplex convergence reached (options.TolFun=' ...
              num2str(options.TolFun) ')' ];
  end
  if istop
    break
  end
end   % while

pars=Simplex('centroid'); % obtain the final value (average).
fval=feval(fun, pars);
funcount=funcount+1;
exitflag=istop;

% output results --------------------------------------------------------------
if istop==0, message='Algorithm terminated normally'; end
output.iterations = iters;
output.algorithm  = [ 'Nelder-Mead Simplex minimization [' mfilename ']' ];
output.message    = message;
output.funcCount  = funcount;

if (istop & strcmp(options.Display,'notify')) | ...
   strcmp(options.Display,'final') | strcmp(options.Display,'iter')
  fmin_private_disp_final(output.algorithm, output.message, output.iterations, ...
    output.funcCount, fun, pars, fval);
end

function [Pars, t]=Simplex(y, StartPars, Steps);
% Nelder-Mead Simplex minimization, implemented as a state machine.
% Usage:
% [Pars, t]=Simplex('init', StartPars, Steps);  % Initialization
% [NewPars, t]=Simplex(y);                      % Iteration
% ok=Simplex('converged', epsilon);             % Test for convergence
% FinalPars=Simplex('centroid');                % Obtain final value
%
% StartPars is the vector of starting parameter values.
% Steps is a vector of initial step sizes (or a scalar which multiplies
% StartPars).
% Pars is the returned array of parameters to try next.
% t is a structure containing the internal state.
% 
% F. Sigworth, 15 March 2003
% Based on a Modula-2 implementation by S. H. Heinemann, 1987 which
% in turn was based on M. Caceci and W. Cacheris, Byte, p. 340, May 1984.
% 
% % Example of use:  minimize (p-q).^6
% p=[2 2 1.5 2]';  % inital guess
% q=[1 1 1 1]';       % true values
% [p,t]=Simplex('init',p);
% iters=0;
% for i=1:200
%     y=sum((p-q).^6);
%     [p,t]=Simplex(y);
%     if ~mod(i, 10)  % every 10 iterations print out the fitted value
%         p'
%     end;
% end;
% p=Simplex('centroid'); % obtain the final value.

persistent PState;
% We make a temporary copy of our state variables for local use.
% This will allow the function to be interrupted without corrupting the
% state variables.
t=PState;   
% The structure elements are the following:
%   t.n  number of elements in the parameter vector
%   t.prow row vector of parameters presently being tested
%   t.simp the simplex matrix, (n+1) row vectors
%   t.vals the (n+1) element column vector of function values
%   t.high index of the highest t.vals element
%   t.low index of the lowest t.vals element
%   t.centr centroid row vector of the best n vertices
%   t.index counter for loops
%   t.state the state variable of the machine

% Interpret the first argument.
if ischar(y)  % an option?
    switch lower(y)
        
        case 'init'  % Initialize the machine, with StartPars being the anchor vertex.
            t.n=prod(size(StartPars));
            t.prow=reshape(StartPars,1,t.n);
            
            % Handle defaults for the step size.
            if nargin <3  % No step size given
                Steps = 0.1;
            end;
            if prod(size(Steps))<t.n  % Only a scalar given.
                zerostep = 1e-3;
                Steps=Steps(1)*t.prow+zerostep*(t.prow==0);
            end;

            % The simplex is (n+1) row vectors.
            t.simp=repmat(t.prow,t.n+1,1);
            for i=1:t.n
                t.simp(i+1,i)=t.simp(i+1,i)+Steps(i);
            end;

            % vals is a column vector of function values
            t.vals=zeros(t.n+1,1);
            
            % Initialize the other variables
            t.index=1;
            t.state=1;
            t.high=0;
            t.low=0;
            t.centr=t.prow;
            Pars=t.prow';
            PState=t;
            
        case 'centroid'  % Return the centroid of the present simplex
            Pars=(sum(t.simp)/(t.n+1))';

        case 'converged'  % Do a convergence test on the vals array.
            err=max(t.vals)-min(t.vals);
            Pars=(t.state==3)&&(err < StartPars);
            
        otherwise
            error('Simplex: unrecognized option');
    end; % switch
    
    
else  % y has a numeric value: this is running mode
    switch t.state
        
        case 1 % Start-up
            t.vals(t.index)=y;  % pick up the function value from last time.
            t.index=t.index+1;
            if t.index <= t.n+1  % continue to fill up the simplex
                t.prow=t.simp(t.index,:);
            else  % Simplex is full, make the first move
                t.state=3;
            end;
            
        case 3  % Test a new vertex
            i=t.high;
            if y < t.vals(i)  % The new vertex is better than some.
                t.simp(i,:)=t.prow;  % replace the worst one.
                t.vals(i)=y;
                if y < t.vals(t.low)  % The new vertex is better than the best,
                   t.prow=t.simp(i,:)+1.1*(t.simp(i,:)-t.centr); % so, expand in the new direction.
                   t.prevy=y;
                   t.state=4;
                else
                    t.state=3;
                end;
            else  % the new vertex is worse than the worst: contract.
                t.prow=0.5*(t.simp(t.high,:)+t.centr);
                t.state=5;
            end;
            
        case 4 % Test an expansion
            if y < t.prevy %t.vals(t.low)  % Accept the expansion
                t.simp(t.high,:)=t.prow;
                t.vals(t.high)=y;
            end;
            t.state=3;
            
        case 5 % Test a contraction
            if y<t.vals(t.high) % Accept the contraction
                t.simp(t.high,:)=t.prow;
                t.vals(t.high)=y;
                t.state=3;
            else %  contract the whole simplex toward the best vertex.
                t.index=1;
                t.simp(1,:)=0.5*(t.simp(1,:)+t.simp(t.low,:));
                prow=t.simp(1,:);
                t.state=6;
            end;
            
        case 6
            t.vals(t.index)=y;  % pick up the function value.
            t.index=t.index+1;
            i=t.index;
            if i <= t.n+1
                t.simp(i,:)=0.5*(t.simp(i,:)+t.simp(t.low,:));
                t.prow=t.simp(i,:);
                t.state=6;  % continue evaluating the vertices
            else
                t.state=3;  % 
            end;
    end; % switch
        
    if t.state==3  % Normal exit mode: sort the vertices and try a reflection.
        % assign min and max
        [x,ind]=sort(t.vals);
        t.low=ind(1);
        t.high=ind(t.n+1);
        % find the excluded centroid
        t.centr=(sum(t.simp)-t.simp(t.high,:))/(t.n);
        % reflect about the centroid from the highest vertex
        t.prow=2*t.centr-t.simp(t.high,:);
    end;
    
    % Copy the output and persistent variables.
    Pars=t.prow';
    PState=t;
end;

