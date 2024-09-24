%==========================================================================
% PROJECT FOR NMoS: INVESTIGATION OF THE INFLUENCE OF AN ADJUSTABLE 
%                   STEP-SIZE IN THE BACKWARD EULER METHOD
%--------------------------------------------------------------------------
%Three methods of adjusting h are compared in their computational time
%accuracy and stability:
%   1. constant h 
%   2. using the second derivative
%   3. adjusting the h in a way to minimize the MSE by using gradient
%       descent
%
%To test these methods functions with a known solution are tested.
%--------------------------------------------------------------------------
%In this file the Implicit Euler MEthod (IEM) is done with second
%derivative
%
%The idea is that if the first derivative is small comapared to the second
%derivative the graph goes near a minima or maxima and therefore 
%adjustment of h is needed to reduce the error and increase stability
%==========================================================================

%--------------------------------------------------------------------------
% PUT IN INITIAL VALUES, FUNCTION AND BOUNDARIES
%--------------------------------------------------------------------------
clear all;
x_start = 5;
t_start = 0;
h_max = 0.1;
t_end = 5;
step_newt = 1000;

%--------------------------------------------------------------------------
% START SOLVING WITH CONSTANT H
%--------------------------------------------------------------------------
t_array(1) = t_start;
i = 1;

while t_array(i)<=t_end %LOOP FOR EULER STEPS
   
   x_sol(1) = x_start;
   f_prime = dx_dt_prime(@dx_dt,t_array(i),x_sol(i));
   f1(i)=f_prime;
   f_2prime = (dx_dt_prime(@dx_dt,t_array(i),x_sol(i)+0.1)-dx_dt_prime(@dx_dt,t_array(i),x_sol(i)))/0.1;
   f2(i)=f_2prime;
   if f_2prime==0
       der2_h=h_max;
   end
   der2_h = abs(f_prime/f_2prime)*0.001;
   if der2_h > h_max || der2_h == 0
       der2_h=h_max;
   end
   der(i)=der2_h;
   x_guess(1)= x_sol(1)+der2_h;%guess of x for the newton method
   x_sol(i+1) = newton(@dx_dt,x_sol(i),x_guess(i),t_array(i)+der2_h,der2_h,step_newt);
   t_array(i+1)=t_array(i)+der2_h;
   x_guess(i+1) = x_sol(i+1)+(x_sol(i+1)-x_sol(i));%guess of x for the newton method
   i = i+1;
end %LOOP FOR EULER STEPS

%--------------------------------------------------------------------------
% REAL VALUES OF X AND MSE
%--------------------------------------------------------------------------

for j=1:i-1 %LOOP FOR REAL VALUES
   x_true(1) = x_start; %solution of x with the newton method 
   mse_value(1) = 0;
   x_true(j+1) = dx_dt_sol(t_array(j+1),x_start);
   mse_value(j+1) = MSE(x_sol(j+1),x_true(j+1));
end %LOOP FOR REAL VALUE
mse_sum = sum(mse_value);
%--------------------------------------------------------------------------
% Visualisation of the Solution
%--------------------------------------------------------------------------
tiledlayout(1,2)
nexttile
plot(t_array,x_sol,"--",'Color',[0, 0.5, 0.3],'LineWidth',1);
hold on
plot(t_array,x_true,"--",'Color',[0.5, 0, 0.3],'LineWidth',1);
xlabel("t");
ylabel("x");
nexttile
plot(t_array,mse_value,"--",'Color',[0.3, 0, 0.5],'LineWidth',1);

%==========================================================================
% FUNCTIONS
%==========================================================================

%--------------------------------------------------------------------------
% FUNCTION TO SOLVE AND FUNCTION FOR REAL SOLUTION (CHANGE THESE FOR YOUR
% PuRPOSE)
%--------------------------------------------------------------------------
function func = dx_dt(t,x)
  func = 2*x;
end

function f_solution = dx_dt_sol(t,x_0)
    f_solution = x_0*exp(2*t);
end

%--------------------------------------------------------------------------
% FUNCTION FOR GETTING DERIVATIVE NUMERICALLY
%--------------------------------------------------------------------------
function func_prime = dx_dt_prime(funct,t,x)
prime_step = 0.0001;
func_prime = (funct(t,x+prime_step)- funct(t,x))+prime_step;

end

%--------------------------------------------------------------------------
% FUNCTION FOR NEWTON METHOD
%--------------------------------------------------------------------------
function x_newt = newton(dx_dt,x_i, x_guess, t, h,step_newt)
  prime_factor = 0.00000001; %For getting the drivative
  error = 0.0001; %Tolerance for finding the root
  %First newton step
  x_search(1) = x_guess;
  
  funct = dx_dt(t,x_guess);
  newton =  x_i + h * funct - x_guess;
  primestep = (x_i + h * dx_dt(t,x_guess+prime_factor) - (x_guess+prime_factor));
  prime = (primestep-(newton))/prime_factor;
  j = 2;
  x_search(j) = x_search(1)-newton/prime;
  while (abs(x_search(j)-x_search(j-1))>error) & j < step_newt %LOOP FOR NEWTON STEPS
    funct = dx_dt(t,x_search(j));
    newton =  x_search(j) + h * funct - x_search(j);
    prime = ((x_i + h * dx_dt(t,x_search(j)+prime_factor) - (x_search(j)+prime_factor))-(newton))/prime_factor;
    x_search(j+1)=x_search(j)-newton/prime;
    j=j+1;
  end%LOOP FOR NEWTON STEPS
  if abs(x_search(j)-x_search(j-1))>error%IF STATEMENT FOR NEWTONSTEPS (if not successful)
      x_newt  = 0;
  end%IF STATEMENT FOR NEWTONSTEPS
  x_newt  = x_search(j); 
end

%--------------------------------------------------------------------------
% FUNCTION FOR MEAN SQUARE ERROR
%--------------------------------------------------------------------------
function MSE = MSE(x_estimate,x_true)
    MSE = (x_estimate-x_true)^2;
end

