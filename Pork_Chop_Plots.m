clear
clc
                        %%%%%%%%%Inputs%%%%%%%%%

z_input = 1.5;                          %Initial guess for value of z
mu_Sun   = 1.32712428e11;               %mu_Sun

                        %%%%%%%%%Pork chop plot%%%%%%%%%

%% Step 1: Choose departure and arrival planets
depPlanet = 'Earth';
arrPlanet = 'Mars';


%% Step 2: Choose an optimal departure and arrival day
day_dep = datetime('20-06-2005','InputFormat','dd-MM-yyyy');
day_arr = datetime('01-12-2008','InputFormat','dd-MM-yyyy');
JD_dep  = juliandate(day_dep);                                  %Convert to Julian Days
JD_arr  = juliandate(day_arr);


%% Step 3: Define the desired time window for arrival and departure (to be added to date above)
tWindowDep = 400; % days added to nominal departure date
tWindowArr = 400; % days added to nominal arrival date


%% Step 4: Define the time steps to build array (defines grid resolution)
tStepDep   = 10;  % departure resolution
tStepArr   = 10;  % arrival resolution


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
[rArray_dep, vArray_dep]=GenerateEphemerides(JDArrayDep, depPlanet);
[rArray_arr, vArray_arr]=GenerateEphemerides(JDArrayArr, arrPlanet);


%% Step 9: Go through all departure/arrival times, and use Lamberts algorithm
% to find needed delta-v's at departure and arrival
fprintf('\n Now building the porkchop plot... ')
counter_3 = 0;

for i = 1:length(JDArrayDep)

    JDi = JDArrayDep(i);

    for j = 1:length(JDArrayArr)

        JDf = JDArrayArr(j);
        
        %Build meshes for contour plot axis
        deltDepMesh(i,j) = JDi-JD_dep;
        deltArrMesh(i,j) = JDf-JD_arr;

        
        % Compute heliocentric orbital velocity at departure and arrival using Lambert's method
        TOF = 86400.0*(JDf - JDi);                                                              % time of flight, in seconds
        [v1Vec,v2Vec] = Function_Lambert_Solver(z_input, mu_Sun, rArray_dep(i,:), rArray_arr(j,:), TOF, 'pro');       
        counter_3 = counter_3 + 1

        % Compute excess velocity for departure and arrival
        vInf_dep = norm(v1Vec - vArray_dep(i,:)); 
        vInf_arr = norm(v2Vec - vArray_arr(j,:));

        % Compute porkchop plot values
        TOFarray(i,j) = JDf - JDi;  
        C3(i,j) = vInf_dep^2;               %C3 = specific energy
        vInf(i,j) = vInf_arr;               %vInf = excess velocity
    end
end
clear i j
fprintf('DONE!\n')

%% PLOT THE PORKCHOP
% 

C3_levels    = linspace(80, 1.8e3, 20);
TOF_levels   = linspace(750, 1700, 20);
V_inf_levels = linspace(13, 33, 20);

col1 = [0.8,0.2,0.2];
col2 = [0.2,0.2,0.8];
col3 = [0.4,0.4,0.4];

% contour sintax: 
% First two arguments are used to define the a grid in the x-y plane
% The third argument is what you want to plot (i.e. Excess velocity)
% The forth argument is the levels in which you want to do them (syntax it
% for better explanation)

close all
figure
set(gcf, 'color', 'w')
hold on
[c1,h1]=contour(deltDepMesh, deltArrMesh, vInf, V_inf_levels, 'showtext','on', 'color', col1,'linewidth',1.5);
[c2,h2]=contour(deltDepMesh, deltArrMesh, C3, C3_levels, 'showtext','on','color', col2,'linewidth',1.5);
[c3,h3]=contour(deltDepMesh, deltArrMesh, TOFarray, TOF_levels, 'showtext','on','color', col3);
hold off
box on
xlabel(['Departure (Days past ', depStr{1},')'],'FontSize',18)
ylabel(['Departure (Days past ', arrStr{1},')'],'FontSize',18)
title([depPlanet, '-to-', arrPlanet, ' Trajectories'],'FontSize',18)
legend({'v_{\infty}','C3','TOF'},'Location','northeastoutside','fontsize',16)



%--------------------------------------------------------------------%
                            %FUNCTIONS%
%% FUNCTION 1: Generating the ephemeris of all of the planets

function [ rArray, vArray ] = GenerateEphemerides( JDArray, planet )
%GenerateEphemerides Generate array of planetary ephemerides

   r_array = zeros(length(JDArray),3);
   v_array = zeros(length(JDArray),3);

   
   for i = 1:length(JDArray)
      [r_array(i,:), v_array(i,:)] = planetEphemeris(JDArray(i), 'Sun', planet,'421');
      % (if you have forgotten what the inputs mean, look at the documentation of 'planetEphemeris')
   end

   rArray = r_array;
   vArray = v_array;

end
