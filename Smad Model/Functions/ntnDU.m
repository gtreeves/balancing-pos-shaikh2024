function [x,J,f] = ntnDU(fHandle,x0,varargin) 
%Newton's method for multivariate functions.
%
%function [x,J] = ntnDU(fHandle,x0,varargin) 
%
%'fHandle' is a function handle that points to the 
%   multivariate function for which the zero is to be found.
%'x0' is an intial guess.
%Optional argument varargin can consist of six things, in this order:
%   (1) 'convergenceCriterion': The convergence tolerance.  Default is 
%       set to 0.001.
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument.  
%	(2) 'nStepsMax': The maximum number of steps before the program quits.
% 		Default is set to 100. 
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument. 
%	(3) 'jacHandle': the Jacobian matrix handle, if the user supplies
%		one.  It is important that the parameters passed to this
%		function are identical to those passed to the fHandle.
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument. 
%	(4)	'dxj': The finite difference size for the numerical Jacobian evaluation in
%		calling jacDU(fHandle,x0,varargin).  Of course, you would leave this out if
%		you specified a 'jacHandle'.  Default is 0.01.
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument.
%	(5)	'alpha': a relaxation constant, that premultiplies the Jacobian inverse.
%		Must be less than one.  Default is 1.
%		If this is not specified, but you still want to specify other 
%		arguments, put an empty string -- '' -- in place of this argument.
%	(6) 'p': A list (which becomes a cell array) of any other parameters that
%		may change the evaluation of the function, but are not actually 
%		varied to find the root.
%       
% x is the final answer, J is the last jacobian, and f is the last function
% evaluation
%warning('off') %raz
nArg = size(varargin,2);
if nArg > 0 && isnumeric(varargin{1})
	convergenceCriterion = varargin{1};
else
	convergenceCriterion = 0.001;
end

if nArg > 1 && isnumeric(varargin{2})
	nStepsMax = varargin{2};
else
	nStepsMax = 100;
end

if nArg > 2 && isnumeric(varargin{3})
	jacHandle = varargin{3};
else
	jacHandle = @jacDU;
end

if nArg > 3 && isnumeric(varargin{4})
	dxj = varargin{4};
else
	dxj = 0.01;
end

if nArg > 4 && isnumeric(varargin{5})
	Alpha1 = varargin{5};
else
	Alpha1 = 1;
end

if nArg > 5
	p = {varargin{6:nArg}};
% 	existOtherParameters = 1;
% else
%     existOtherParameters = 0;
end

% if existOtherParameters
    %This part of the if statement is executed if other parameters are passed to the function fHandle.
	
% 	f = fHandle(x0,p{:});
%     J = jacHandle(fHandle,x0,f,dxj,p{:}); %Initializing the code.
% 	
%     x = x0 - Alpha1*(J \ f);
%     delta_x = norm(x - x0);
%     x0 = x;
%     nSteps = 1;

	f = fHandle(x0,p{:});
	J = jacHandle(fHandle,x0,f,dxj,p{:});
	delta_x = 1;
	nSteps = 0;

    %The MAIN LOOP.
    while (norm(f) > Alpha1*convergenceCriterion) && (nSteps < nStepsMax) && (rcond(J) > 0.01)
		f = fHandle(x0,p{:});
        J = jacHandle(fHandle,x0,f,dxj,p{:});
		delta_x = - Alpha1*(J \ f);
        x = x0 + delta_x;
%         delta_x = norm(x - x0);
        x0 = x;
        nSteps = nSteps + 1;
    end
    
    x = x0;
% else
%     %This part of the if statement (the 'else' part) is executed if no other parameters are passed to the function fHandle.
%     f = feval(fHandle,x0);
% 	J = feval(jacHandle,fHandle,x0,f); %Initializing the code.
%     x = x0 - Alpha1*(J \ f);
%     delta_x = norm(x - x0);
%     x0 = x;
%     nSteps = 1;
% 
%     %The MAIN LOOP.
%     while (delta_x > Alpha1*convergenceCriterion) & (nSteps < nStepsMax)
%        f = feval(fHandle,x0);
% 	   J = feval(jacHandle,fHandle,x0,f,dxj);
%         x = x0 - Alpha1*(J \ f);
%         delta_x = norm(x - x0);
%         x0 = x;
%         nSteps = nSteps + 1;
%     end    
% end
%if delta_x > Alpha1*convergenceCriterion
%    x = 'Did not converge';
end