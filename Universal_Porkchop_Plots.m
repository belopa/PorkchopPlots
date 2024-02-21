clear
clc

%Function used: Function_Lambert_Solver
                        %%%%%%%%%Inputs%%%%%%%%%

z_input       = -1.6e3;                   %Initial guess for value of z

%Planetary constants
mu_Sun        = 1.32712428e11;            %mu_Sun [km^3/s^2]
mu_Mercury_km = 2.203e4;                  %Standard gravitational parameter of Mercury [km^3/s^2]
mu_Venus_km   = 3.249e5;                  %Standard gravitational parameter of Venus [km^3/s^2]
mu_Earth_km   = 3.986e5;                  %Standard gravitational parameter of Earth [km^3/s^2]
mu_Mars_km    = 4.283e4;                  %Standard gravitational parameter of Mars [km^3/s^2]
mu_Jupiter_km = 1.267e8;                  %Standard gravitational parameter of Jupiter [km^3/s^2]
mu_Saturn_km  = 3.794e7;                  %Standard gravitational parameter of Saturn [km^3/s^2]
mu_Uranus_km  = 5.795e6;                  %Standard gravitational parameter of Uranus [km^3/s^2]
mu_Neptune_km = 6.837e6;                  %Standard gravitational parameter of Neptune [km^3/s^2]
d_SunEarth    = 150e6;                    %Distance between the Sun and Earth [km]
d_SunMars     = 228e6;                    %Distance between the Sun and Mars [km]
r_Earth       = 6378.1;                   %Radius of Earth [km]
r_Mars        = 3389.5;                   %Radius of Mars [km] 

                        %%%%%%%%%Pork chop plot%%%%%%%%%

%% Step 1: Choose departure and arrival planets and their properties
%For convenience, create new variable for the departure and arrival planets
%(so you don't have to change every variable every time you change the planet of
%interest)
depPlanet = 'Earth';
arrPlanet = 'Mars';
dep_mu    = mu_Earth_km;                %Standard gravitational parameter of departure planet [km^3/s^2]          
arr_mu    = mu_Mars_km;                 %Standard gravitational parameter of arrival planet [km^3/s^2]
d_SunP1   = d_SunEarth;                 %Distance of departure planet to Sun [km]
d_SunP2   = d_SunMars;                  %Distance of arrival planet to Sun [km]
r_P1      = r_Earth;                    %Radius of departure planet [km]
r_P2      = r_Mars;                     %Radius of arrival planet [km]
alt_P1    = 300;                        %Altitude of initial orbit around departure planet [km]
alt_P2    = 200;                        %Altitude of final orbit around departure planet [km]



%% Step 2: Choose an optimal departure and arrival day
day_dep = datetime('04-08-2026','InputFormat','dd-MM-yyyy');
day_arr = datetime('20-03-2027','InputFormat','dd-MM-yyyy');
%Convert to Julian Days
JD_dep  = juliandate(day_dep);                                  
JD_arr  = juliandate(day_arr);
%Define the desired time window for arrival and departure
tWindowDep = 180;                                      
tWindowArr = 360;                                     
%Define the time steps to build array (defines grid resolution)
tStepDep   = 1;                                         
tStepArr   = 2;                                        
%Print what the arrival and departure days are (useful when plotting later)
depDate=datetime(JD_dep,'convertfrom','juliandate','Format','dd-MMM-yyy');
depStr=cellstr(depDate);
arrDate=datetime(JD_arr,'convertfrom','juliandate','Format','dd-MMM-yyy');
arrStr=cellstr(arrDate);
%Create an array of departure and arrival windows in Julian Days
JDArrayDep = JD_dep : tStepDep : JD_dep+tWindowDep;
JDArrayArr = JD_arr : tStepArr : JD_arr+tWindowArr;
%Make sure departure and arrival windows dont overlap
if JDArrayDep(end) >= JDArrayArr(1)
    fprintf(' Error: Some arrival times before departure times! \n');
    fprintf('        Try adjusting arrival date and tWindow \n');
    return
end

%% Step 3: Generate array of planetary ephemerides at each departure/arrival time
% (The function used (GenerateEphemerides) can be found at the bottom of
% the doc)
[rArray_dep_km, vArray_dep_km]=GenerateEphemerides(JDArrayDep, depPlanet);
[rArray_arr_km, vArray_arr_km]=GenerateEphemerides(JDArrayArr, arrPlanet);


%% Step 4: Create cells to have all three velocity components in the same cell for correct indexing below
%Empty cells
cell_arr_planet = cell(length(JDArrayDep),length(JDArrayDep));
cell_dep_planet = cell(length(JDArrayDep),length(JDArrayDep));

%Input velocities into empty cells:
%Want to create a departure array with the same velocity in each COLUMN for
%the correct indexing in loop below (Step 5)
for i = 1:length(JDArrayDep)
    for j = 1:length(JDArrayDep)
        cell_dep_planet{j,i} = vArray_dep_km(i,:);
    end
end

clear i j

%Want to create an arrival array with the same velocity in each ROW for
%the correct indexing in loop below (Step 5)
for j = 1:length(JDArrayArr)
    for i = 1:length(JDArrayArr)
        cell_arr_planet{j,i} = vArray_arr_km(j,:);
    end
end

clear i j

%Rename
Velocity_Planet1 = cell_dep_planet;
Velocity_Planet2 = cell_arr_planet;

%% Step 5: Go through all departure/arrival times, and use Lamberts algorithm
% to find needed delta-v's at departure and arrival
fprintf('\n Now building the porkchop plot... ')
counter_3 = 0;

for i = 1:length(JDArrayDep)

    for j = 1:length(JDArrayArr)
          
        % Compute heliocentric orbital velocity at departure and arrival using Lambert's method
        TOF(i,j)                = 86400.0*(JDArrayArr(:,j) - JDArrayDep(:,i));         % time of flight, in seconds 
        [v1Vec{i,j},v2Vec{i,j}, Elements] = Function_Lambert_Solver(z_input, mu_Sun, rArray_dep_km(i,:), rArray_arr_km(j,:), TOF(i,j), 'pro');       
        counter_3               = counter_3 + 1

    end
end

clear i j

% Intent to get rid of imaginary values:
% M0 = Elements(6);
% mask = imag(M0) ~= 0;
% v1Vec_real = v1Vec(~mask, :);

%% Calculate excess velocity & specific energy & time of flight
%Create cells for correct indexing
cell_vinf1 = cell(length(JDArrayDep),length(JDArrayDep));
cell_vinf2 = cell(length(JDArrayDep),length(JDArrayDep));

%Rename
vInf_dep = cell_vinf1;
vInf_arr = cell_vinf2;

%Compute excess velocity at arrival and departure using correct indexing
for i = 1:length(JDArrayDep)
    for j = 1:length(JDArrayDep)
        vInf_dep{i,j} = v1Vec{i,j} - Velocity_Planet1{i,j}; 
        vInf_dep_array(i,j) = norm(v1Vec{i,j} - Velocity_Planet1{i,j});
        vInf_arr{i,j} = v2Vec{i,j} - Velocity_Planet2{i,j};
        vInf_arr_array(i,j) = norm(v2Vec{i,j} - Velocity_Planet2{i,j});
    end
end

clear i j

%Rename for plotting
TOF_Array     = TOF./86400;                   %Time of flight in days
C3            = vInf_dep_array.^2;            %Specific energy
vInf          = vInf_arr_array;               %Excess velocity


%Calculating total delta V assuming Hohmann transfer 
v_apoapsis_hyperbola     = sqrt(vInf_dep_array.^2+((2*dep_mu)/alt_P1));        %Velocity of hyperbolic orbit closest to Planet 1
v_circ_wrtP1             = sqrt(dep_mu/alt_P1);                                %Velocity of circular orbit around P1
DELTAV_1                 = v_apoapsis_hyperbola-v_circ_wrtP1;                  %Delta V required for transfer

a_transfer               = -arr_mu./(vInf.^2);                                 %Semi-major axis of hyperbolic transfer orbit
v_periapsis_hyperbola    = sqrt(((2*arr_mu)/alt_P2)-(arr_mu./a_transfer));      %Velocity of hyperbolic orbit closest to Planet 2
v_circ_wrtP2             = sqrt(arr_mu/alt_P2);                                %Velocity of circular orbit around P2
DELTAV_2                 = v_periapsis_hyperbola-v_circ_wrtP2;                 %Delta V required for transfer

TOTALDELTAV              = abs(DELTAV_1)+abs(DELTAV_2);                        %Total delta V


%% PLOTTING POSITION OF PLANETS

figure(1)
grid on
plot3(rArray_arr_km(:,1),rArray_arr_km(:,2),rArray_arr_km(:,3),'o')
title(['Orbit of ' depPlanet ' and ' arrPlanet ' around Sun'])
hold on
plot3(rArray_dep_km(:,1),rArray_dep_km(:,2),rArray_dep_km(:,3),'o')
legend([arrPlanet ' orbit'],[depPlanet ' orbit'])
grid off
hold off

fprintf('Position plot: DONE!\n')

%% PLOT THE PORKCHOP

figure(2)
%hold on
contour(JDArrayDep,JDArrayArr, vInf, 20,'ShowText','on','color','b');
%contour(JDArrayDep, JDArrayArr, C3, 20, 'color','r','ShowText','on');
%contour(JDArrayDep, JDArrayArr, TOF_Array, 10, 'color', 'g','ShowText','on', 'LineStyle','--', 'LineWidth',1.5);
box on
%hold off
xlabel(['Departure (Days past ', depStr{1},')'],'FontSize',18)
ylabel(['Arrival (Days past ', arrStr{1},')'],'FontSize',18)
title([depPlanet, '-to-', arrPlanet, ' Trajectories'],'FontSize',18)
%legend({'v_{\infty}','C3','TOF'},'Location','northeastoutside','fontsize',16)

fprintf('Porkchop plot: DONE!\n')


% Plotting total delta V 

figure(3)
grid on
contourf(JDArrayDep,JDArrayArr, vInf, 40,'ShowText','on');
box on
xlabel(['Departure (Days past ', depStr{1},')'],'FontSize',18)
ylabel(['Arrival (Days past ', arrStr{1},')'],'FontSize',18)
title(['vInf: depPlanet:', depPlanet, '-to-', arrPlanet, ' Trajectories'],'FontSize',18)

fprintf('Porkchop plot of delta V: DONE!\n')

%% SURFACE PLOT
figure(4)
surf(JDArrayDep,JDArrayArr,C3)
xlabel(['Departure (Days past ', depStr{1},')'],'FontSize',18)
ylabel(['Arrival (Days past ', arrStr{1},')'],'FontSize',18)
zlabel('v_{\infty}','FontSize',18)
title([depPlanet, '-to-', arrPlanet, ' Trajectories'],'FontSize',18)

fprintf('Surface plot: DONE!\n')
%--------------------------------------------------------------------%
                            %FUNCTIONS%
%% FUNCTION 1: Generating the ephemeris of all of the planets

function [ rArray, vArray ] = GenerateEphemerides( JDArray,  planet)

   r_array = zeros(length(JDArray),3);
   v_array = zeros(length(JDArray),3);

   
   for i = 1:length(JDArray)
      [r_array(i,:), v_array(i,:)] = planetEphemeris(JDArray(i),'Sun',planet,'421');
      % (if you have forgotten what the inputs mean, look at the documentation of 'planetEphemeris')
   end

   rArray = r_array;
   vArray = v_array;

end
