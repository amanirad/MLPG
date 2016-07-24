%
% What should we give in? An initial condition in terms of something.
% An f-function and a g-function whatever they are. All data must be input.
% Then S of course.

%
% Check up sending function handles as arguments
%
function [lambda,mu]=BDF2LSSolver(S,k,beta1,beta2,u0e,BCrhs)
%
%
  %  
  % Initialize all previous coefficients to zero
  %
  lambda_1 = 0; lambda_2 = 0; mu_1 = 0; mu_2 = 0;
  
  M = length(k); % The number of time steps 
  for step = 1:M
    %
    % If time does not start at zero, the offset must be incorporated
    % into the functions. Functions should always have the form AArhs(x,t)
    %
    time(step) = sum(k(1:step));
    %
    % Compute the right hand side for the current step
    %
    g = BCrhs(time(step));
    %
    % Compute the BDF2-right hand side beta_1*u_{n-1}-beta_2*u_{n-2}
    % The first two steps are special cases.
    %
    if (step==1)
      f = u0e;
    elseif (step==2)
      f = beta1(step)*(S.Ael*lambda_1+S.Aem*mu_1)-beta2(step)*u0e;
    else
      f = S.Ael*(beta1(step)*lambda_1-beta2(step)*lambda_2) + ...
          S.Aem*(beta1(step)*mu_1    -beta2(step)*mu_2);
    end
    %
    % Solve step 1: Abm*w=g^n
    %
    w = S.Ub\(S.Lb\(S.Pb*g));
    %
    % Solve step 2: S*lambda^n = f^n-Cem*w
    %
    f = f - S.Cem*w;
    b = S.Q'*f;
    lambda = S.R\b;
    %
    % Solve step 3: Abm*v = Abl*lambda^n
    %
    v = S.Ub\(S.Lb\(S.Pb*(S.Abl*lambda)));
    %
    % Solve step 4: mu^n = w-v;
    %
    mu = w-v;
    %
    % Move the values from the previous step
    %
    lambda_2 = lambda_1;
    mu_2 = mu_1;
    lambda_1 = lambda;
    mu_1 = mu;
    %
    % If we want to compute the mass, we need location and ep.
    %
    %coeff_l = 1 - normcdf(0, xc,  1/sqrt(2)/ep(1));
    %coeff_m = 1 - normcdf(0, xc,  1/sqrt(2)/ep(1));
    %mass(step)=sqrt(pi)/ep(1)*(sum(coeff_l.*lambda{1})+sum(coeff_m.*mu{1}));
  end
