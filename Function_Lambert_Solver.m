function [Velocity_Vector1,Velocity_Vector2, Elements] = Function_Lambert_Solver(z, mu,R1,R2,R2_inv,delta_t,string)

                        %%%%%%%%%INTRODUCTION%%%%%%%%
%%In this function we want to calculate the velocity vectors of a point on
%%a path from point 1 (with position R1) to point 2 (with position R2). The
%%variable delta_t describes the time of flight and the 'string' is used to
%%input whether you want the orbit to be prograde or retrograde.

%%Note that all the functions used throughout are found at the end of this
%%document.

    %%%%%%%%%ALL OF THESE EQUATIONS CAN BE FOUND IN "ORBITAL MECHANICS FOR
                        %%%%%%%%%ENGINEERING STUDENTS" 

%Inputs of funciton:
    % z = initial guess for Newton's iteration
    % mu = Standard gravitational parameter of the Sun [km^3/s^2]
    % R1 = position of departure planet
    % R2 = position of arrival planet 
    % delta_t = tof (time of flight)
    % string = 'pro' or 'retro'
%Outputs of function:
    % Velocity_Vector1 = velocity required by satellite at departure planet
    % Velocity_Vector2 = velocity required by satellite at arrival planet

%% Step 1:  Calculate the norm, dot and cross product of the position vectors
r1 = norm(R1);
r2 = norm(R2);
dot_product = dot(R1,R2);
cross_product = cross(R1,R2);

%% Step 2: Calculate the change in true anomaly (sigma)
% Checking whether the orbit is prograde or retrograde to change the TA accordingly (Eq
% 5.26)

if strcmp(string, 'pro')

    if cross_product(:,3) >= 0
        sigma =  acos((dot_product)./(r1*r2));
    else  
        sigma = (2*pi) - acos((dot_product)./(r1*r2));
    end

elseif strcmp(string, 'retro')

    if cross_product(:,3) < 0
        sigma =  acos((dot_product)./(r1*r2));
    else  
        sigma = (2*pi) - acos((dot_product)./(r1*r2));
    end

end

%% Step 3: Simplification expression  
A = sin(sigma).*sqrt((r1*r2)./(1-cos(sigma)));        % (Eq 5.35)

%% Step 4: Find the rough value of z where F = 0
F_something = F(z,delta_t);                           % Initial value of F using initial z and tof values 

while F_something < 0
    z = z + 0.1;
    F_something = F(z,delta_t);
end

%% Step 5: Calculate what the value of z is using initial guess given in Step 4
% (Newton's procedure Pg 209)

tol       = 1e-3;                    %Tolerance
nmax      = 1000000;                 %Limit on iterations
F_ratio   = F_something/F_der(z);    %Initial value of F_ratio
counter_2 = 0;

while abs(F_ratio) > tol && counter_2 <= nmax
    z         = z - F_ratio;
    F_ratio   = F(z,delta_t)/F_der(z);
    counter_2 = counter_2 + 1;
end

% Display error if the maximum number of iterations is exceeded:
if counter_2 >= nmax
    fprintf('\n\n **Number of iterations exceeds %g in ''lambert'' \n\n ',nmax)
end


%% Step 5: Calculating Lagrange coefficients (Expressed in terms of universal anomaly for this analysis)
f     = 1 - (Y(z)./r1);            % Lagrange function (Eq 5.46a)
g     = A.*sqrt(Y(z)./mu);         % Lagrange function (Eq 5.46b)
g_dot = 1 - (Y(z)./r2);            % Lagrange function (Eq 5.46d)

%% Step 6: Calculating Velocity Vectors (Rearranging the relation between Lagrange coefficients and position&velocity vectors)

Velocity_Vector1 = (1./g).*(R2-(f.*R1));              % Velocity Vector 1 (Eq 5.28)
Velocity_Vector2 = (1./g).*((g_dot.*R2)-R1);          % Velocity Vector 2 (Eq 5.29)

R = [R1 R2];                                        % Position Vector
V = [Velocity_Vector1 Velocity_Vector2];            % Velocity Vector


%% Step 7: Converting Position Vectors into Keplerian Elements 
% https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf

r1     = norm(R1);
v1     = norm(Velocity_Vector1);
h      = cross(R1,Velocity_Vector1);                             % Orbital momentum vector
e_vec  = ((cross(Velocity_Vector1,h))/mu) - (R1/norm(R1));       % Eccentric vector 
e      = norm(e_vec);                                            % e scalar
E      = 2 * atan( (tan(v/2)) / (sqrt((1+e)/(1-e))) );           % Eccentric anomaly
i      = (acos(h(3) / norm(h)));                                 % Inclination
M0     = E - e*sin(E);                                           % Mean anomlay
a      = 1/ ( (2/r1) - (v1^2/mu) );                              % Semi-major axis
T      = 2 * pi * sqrt((a^3)/mu);                                % Period

k     = [0 0 1];
n     = cross(k', h);                                            % Vector that points towards the ascending node

%True anomaly
if dot(R1, Velocity_Vector1) >= 0                              
    v = acos((dot(R1,Velocity_Vector1))/(r1*v1));
else
    v = 2*pi - acos((dot(R1,Velocity_Vector1))/(r1*v1));
end

% Ascending node
if n(2)>= 0 
    O = acos(n(1)/norm(n));                                
else
    O = 2*pi - acos(n(1)/norm(n));
end

%Argument of pariapsis 
if e_vec(3) >= 0 
    w = acos((dot(n,e_vec))/(norm(n)*e));                   
else
    w = 2*pi - acos((dot(n,e_vec))/(norm(n)*e));
end

Elements = [a e i O w M0];


% Intent to get rid of values that give imaginary semi major axis:
% mask = imag(T) == 0;
% T_real = T(mask,:);


%--------------------------------------------------------------------%
                            %FUNCTIONS%

%% FUNCTION 1: Value of y (Equation 5.38)
function Value_y = Y(z)
Value_y = r1 + r2 + (A.*(((z.*S(z))-1)./sqrt(C(z))));
end

%% FUNCTION 2: Value of F(z) function (Equation 5.40)
function F_function = F(z, delta_t)
    F_function = (((Y(z)./C(z)).^1.5).*S(z)) + A.*sqrt(Y(z)) - (sqrt(mu)*delta_t);
end

%% FUNCTION 3: Value of derivative of F(z) function (Equation 5.43)
function Derivative_F = F_der(z)
if z == 0 
    Derivative_F = ((sqrt(2)/40)*(Y(0)^1.5)) + (A/8).*(sqrt(Y(0)) + A.*sqrt(1/(2.*Y(0))));
else
    Derivative_F = ((Y(z)./C(z)).^1.5) * ((1./(2.*z)).*(C(z) - (3.*S(z))/(2.*C(z))) ...
          + (3.*(S(z).^2))./(4.*C(z))) + A/.8.*(((3.*S(z))./C(z)).*sqrt(Y(z)) ...
          + A.*sqrt(C(z)./Y(z)));
end

end

%% FUNCTION 4a: Stumpff function S(z) (Equation 3.49)
function StumpffFunction_S = S(z)
if z>0
    StumpffFunction_S = (sqrt(z)-sin(sqrt(z)))/((sqrt(z))^3);
elseif z < 0 
    StumpffFunction_S = (sinh(sqrt(-z))-sqrt(-z))/((sqrt(-z))^3);
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
