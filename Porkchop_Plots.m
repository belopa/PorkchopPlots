clear
clc
                        %%%%%%%%%Inputs%%%%%%%%%

z_input       = 2;                        %Initial guess for value of z
mu_Sun        = 1.32712428e11;            %mu_Sun (km)
mu_Mercury_km = 2.203e4;                  %Standard gravitational parameter of Mercury
mu_Venus_km   = 3.249e5;                  %Standard gravitational parameter of Venus
mu_Earth_km   = 3.986e5;                  %Standard gravitational parameter of Earth
mu_Mars_km    = 4.283e4;                  %Standard gravitational parameter of Mars
mu_Jupiter_km = 1.267e8;                  %Standard gravitational parameter of Jupiter
mu_Saturn_km  = 3.794e7;                  %Standard gravitational parameter of Saturn
mu_Uranus_km  = 5.795e6;                  %Standard gravitational parameter of Uranus
mu_Neptune_km = 6.837e6;                  %Standard gravitational parameter of Neptune

                        %%%%%%%%%Pork chop plot%%%%%%%%%

%% Step 1: Choose departure and arrival planets
depPlanet = 'Earth';
arrPlanet = 'Mars';
dep_mu    = mu_Earth_km;
arr_mu    = mu_Mars_km;

%% Step 2: Choose an optimal departure and arrival day
day_dep = datetime('01-05-2020','InputFormat','dd-MM-yyyy');
day_arr = datetime('01-12-2020','InputFormat','dd-MM-yyyy');

%Convert to Julian Days
JD_dep  = juliandate(day_dep);                                  
JD_arr  = juliandate(day_arr);


%% Step 3: Define the desired time window for arrival and departure (to be added to date above)
tWindowDep = 200; % days added to nominal departure date
tWindowArr = 200; % days added to nominal arrival date


%% Step 4: Define the time steps to build array (defines grid resolution)
tStepDep   = 1;  % departure resolution
tStepArr   = 1;  % arrival resolution


%% Step 5: Print what the arrival and departure days are in Command Window
% (user interface friendly, not essential)
fprintf('\nGenerating the porkchop plot for trajectories \n')
fprintf('departing %s and arriving at %s.\n',depPlanet,arrPlanet)
depDate=datetime(JD_dep,'convertfrom','juliandate','Format','dd-MMM-yyy');
depStr=cellstr(depDate);
arrDate=datetime(JD_arr,'convertfrom','juliandate','Format','dd-MMM-yyy');
arrStr=cellstr(arrDate);
fprintf('Nominal departure date: %s \n',depStr{1})
fprintf('Nominal arrival date: %s \n \n', arrStr{1})


%% Step 6: Create an array of departure and arrival windows in Julian Days
JDArrayDep = JD_dep : tStepDep : JD_dep+tWindowDep;
JDArrayArr = JD_arr : tStepArr : JD_arr+tWindowArr;


%% Step 7: Test to make sure all arrival times are after departure times
if JDArrayDep(end) >= JDArrayArr(1)
    fprintf('\n*********************************************** \n');
    fprintf(' Error: Some arrival times before departure times! \n');
    fprintf('        Try adjusting arrival date and tWindow \n');
    fprintf('*********************************************** \n');
    return
end


%% Step 8: Generate array of planetary ephemerides at each departure/arrival time
% (The function used (GenerateEphemerides) can be found at the bottom of
% the doc)
[rArray_dep_km, vArray_dep_km]=GenerateEphemerides(JDArrayDep, depPlanet);
[rArray_arr_km, vArray_arr_km]=GenerateEphemerides(JDArrayArr, arrPlanet);

%used later
rArray_arr_km_inv = rArray_arr_km';

%% Create cells to have all three velocity components in the same cell for correct indexing below
%Empty cells
cell_arr_planet = cell(length(JDArrayDep),length(JDArrayDep));
cell_dep_planet = cell(length(JDArrayDep),length(JDArrayDep));

%Input velocities into empty cells
for i = 1:length(JDArrayDep)
    for j = 1:length(JDArrayDep)
        cell_dep_planet{j,i} = vArray_dep_km(i,:);

    end
end

clear i j

for j = 1:length(JDArrayArr)
    for i = 1:length(JDArrayArr)
        cell_arr_planet{j,i} = vArray_arr_km(j,:);

    end
end

clear i j

%Rename
Velocity_Planet1 = cell_dep_planet;
Velocity_Planet2 = cell_arr_planet;



%% Step 9: Go through all departure/arrival times, and use Lamberts algorithm
% to find needed delta-v's at departure and arrival
fprintf('\n Now building the porkchop plot... ')
counter_3 = 0;

for i = 1:length(JDArrayDep)

    JDi = JDArrayDep(i);

    for j = 1:length(JDArrayArr)

        JDf = JDArrayArr(j);
        
        % Build meshes for contour plot axis
        deltDepMesh(i,:) = JDArrayDep(:,i)-JD_dep;
        deltArrMesh(j,:) = JDArrayArr(:,j)-JD_arr;
        
        
        % Compute heliocentric orbital velocity at departure and arrival using Lambert's method
        TOF(i,j)                = 86400.0*(JDArrayArr(:,j) - JDArrayDep(:,i));         % time of flight, in seconds (conversion from days to s)
        [v1Vec{i,j},v2Vec{i,j}, Elements] = Function_Lambert_Solver(z_input, mu_Sun, rArray_dep_km(i,:), rArray_arr_km(j,:), rArray_arr_km_inv(:,j), TOF(i,j), 'pro');       
        counter_3               = counter_3 + 1

      
    end
end

clear i j

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


Array_Dep = datetime(JDArrayDep, 'ConvertFrom', 'juliandate');
Array_Arr = datetime(JDArrayArr, 'ConvertFrom', 'juliandate');



%% PLOTTING POSITION OF PLANETS

figure(1)
plot3(rArray_arr_km(:,1),rArray_arr_km(:,2),rArray_arr_km(:,3),'o')
title(['Orbit of ' depPlanet ' and ' arrPlanet ' around Sun'])
grid on
hold on
plot3(rArray_dep_km(:,1),rArray_dep_km(:,2),rArray_dep_km(:,3),'o')
legend([arrPlanet ' orbit'],[depPlanet ' orbit'])
hold off

fprintf('Position plot: DONE!\n')

%% PLOT THE PORKCHOP

figure(2)
hold on
grid on
contour(JDArrayDep,JDArrayArr, vInf, 20,'ShowText','on','color','r');
contour(JDArrayDep, JDArrayArr, C3, 10, 'color','b','LineWidth',1.5,'ShowText','on');
contour(JDArrayDep, JDArrayArr, TOF_Array, 10, 'color', 'g','ShowText','on', 'LineStyle','--', 'LineWidth',1.5);
box on
hold off
xlabel(['Departure (Days past ', depStr{1},')'],'FontSize',18)
ylabel(['Arrival (Days past ', arrStr{1},')'],'FontSize',18)
title([depPlanet, '-to-', arrPlanet, ' Trajectories'],'FontSize',18)
legend({'v_{\infty}','C3','TOF'},'Location','northeastoutside','fontsize',16)

fprintf('Porkchop plot: DONE!\n')


%% SURFACE PLOT
figure(3)
surf(JDArrayDep,JDArrayArr,vInf)
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
