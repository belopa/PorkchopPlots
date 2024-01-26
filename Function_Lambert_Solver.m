function [Velocity_Vector1,Velocity_Vector2, Elements] = Function_Lambert_Solver(z, mu,R1,R2,R2_inv,delta_t,string)

                        %%%%%%%%%INTRODUCTION%%%%%%%%
%%In this function we want to calculate the velocity vectors of a point on
%%a path from point 1 (with position R1) to point 2 (with position R2). The
%%variable delta_t describes the time of flight and the 'string' is used to
%%input whether you want the orbit to be prograde or retrograde.

%%Note that all the functions used throughout are found at the end of this
%%document. In such way, the inputs (like z or delta_t) can be inputted
%%into the function and used in certain areas of the code without having to
%%include it every time (if that makes sense). The code will refer to the
%%functions included at the bottom every time it needs to and input the
%%according value of z or delta_t that you call.

    %%%%%%%%%ALL OF THESE EQUATIONS CAN BE FOUND IN "ORBITAL MECHANICS FOR
                        %%%%%%%%%ENGINEERING STUDENTS"            

%% Step 1&2:  Calculate the norm of the position vectors to give the "magnitude"
%% Calculate the change in true anomaly 

r1 = norm(R1);
r2 = norm(R2);
dot_product = dot(R1,R2_inv);
multiple_position = norm(R1)*norm(R2_inv);
c12 = cross(R1, R2);
division_position = (dot_product./multiple_position');
sigma = acos(division_position);                         % Change in true anomaly (Eq 5.23) 



% Checking whether the orbit is prograde or retrograde to change the TA accordingly (Eq
% 5.26)
if strcmp(string,'pro')
   if c12(3) <= 0
        sigma = 2*pi - sigma;
   end
elseif strcmp(string,'retro')
   if c12(3) >= 0
        sigma = 2*pi - sigma;
   end
end

%% Step 3: Simplification expression  
A = sin(sigma).*sqrt((r1*r2)./(1-cos(sigma)));        % (Eq 5.35)


%% Step 4: Calculate what the value of z is using initial guess given in Step 4
% (Newton's procedure Pg 209)

tol       = 1e-3;          %Tolerance
nmax      = 100000;        %Limit on iterations
F_ratio   = 0.001;
counter_2 = 0;

while abs(F_ratio) > tol && counter_2 <= nmax
    counter_2 = counter_2 + 1
    F_ratio   = F(z,delta_t)/F_der(z);
    z         = z - F_ratio;
end

% Display error if the maximum number of iterations is exceeded:
if counter_2 >= nmax
    fprintf('\n\n **Number of iterations exceeds %g in ''lambert'' \n\n ',nmax)
end


%% Step 5: Calculating Lagrange function (Expressed in terms of universal anomaly for this analysis)
f     = 1 - (y(z)./r1);            % Lagrange function (Eq 5.46a)
g     = A.*sqrt(y(z)./mu);         % Lagrange function (Eq 5.46b)
g_dot = 1 - (y(z)./r2);            % Lagrange function (Eq 5.46d)

%% Step 6: Calculating Velocity Vectors (Rearranging the relation between Lagrange coefficients and position&velocity vectors)

Velocity_Vector1 = (1./g).*(R2-(f.*R1));              % Velocity Vector 1 (Eq 5.28)
Velocity_Vector2 = (1./g).*((g_dot.*R2)-R1);          % Velocity Vector 2 (Eq 5.29)

R = [R1 R2];                                        % Position Vector
V = [Velocity_Vector1 Velocity_Vector2];            % Velocity Vector


%% Step 7: Converting Position Vectors into Keplerian Elements

r     = norm(R1);
v     = norm(Velocity_Vector1);
v_r   = dot(R1/r, Velocity_Vector1);
v_p   = sqrt(v^2 - v_r^2);
h_vec = cross(R1, Velocity_Vector1);
h     = norm(h_vec);                                                % Angular momentum
i     = rad2deg(acos(h_vec(3) / h));                                % Inclination
K     = [0, 0, 1];  
N_vec = cross(K, h_vec);        
N     = norm(N_vec);    
Omega = 360 - rad2deg(2 * pi - acos(N_vec(1)/N));                   % RAAN
e_vec = (cross(Velocity_Vector1, h_vec) / mu) - (R1 / r);
e     = norm(e_vec);                                                % Eccentricity
w     = 360 - rad2deg(2 * pi - acos(dot(N_vec, e_vec) / (N * e)));  % Argument of Perigee
nu    = 360 - rad2deg(acos(dot(e_vec, R1) / (e * r)));              % Mean Anomaly
a     = h^2/(mu*(1 - e^2));                                         % Semi-major axis

Elements = [a e i Omega w nu];



%--------------------------------------------------------------------%
                            %FUNCTIONS%

%% FUNCTION 1: Value of y (Equation 5.38)
function Value_y = y(z)
Value_y = r1 + r2 + (A.*(((z.*S(z))-1)./sqrt(C(z))));
end

%% FUNCTION 2: Value of F(z) function (Equation 5.40)
function F_function = F(z, delta_t)
    F_function = (((y(z)./C(z)).^1.5).*S(z)) + A.*sqrt(y(z)) - (sqrt(mu)*delta_t);
end

%% FUNCTION 3: Value of derivative of F(z) function (Equation 5.43)
function Derivative_F = F_der(z)
if z == 0 
    Derivative_F = ((sqrt(2)/40)*(y(0)^1.5)) + (A/8).*(sqrt(y(0)) + A.*sqrt(1/(2.*y(0))));
else
    Derivative_F = ((y(z)./C(z)).^1.5) * ((1./(2.*z)).*(C(z) - (3.*S(z))/(2.*C(z))) ...
          + (3.*(S(z).^2))./(4.*C(z))) + A/.8.*(((3.*S(z))./C(z)).*sqrt(y(z)) ...
          + A.*sqrt(C(z)./y(z)));
end

end

%% FUNCTION 4a: Stumpff function S(z) (Equation 3.49)
function StumpffFunction_S = S(z)
if z>0
    StumpffFunction_S = (sqrt(z)-sin(sqrt(z)))/((sqrt(z))^3);
elseif z < 0 
    StumpffFunction_S = (sinh(sqrt(-z))-sqrt(-z))/((sqrt(z))^3);
else
    StumpffFunction_S = 1/6;
end
end

%% FUNCTION 4: Stumpff function S(z) (Equation 3.50)
function StumpffFunction_C = C(z)
if z>0
    StumpffFunction_C = (1-cos(sqrt(z)))/z;
elseif z < 0 
    StumpffFunction_C = (cosh(sqrt(-z))-1)/-z;
else
    StumpffFunction_C = 1/2;
end
end

end

end
