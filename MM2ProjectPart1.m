function internals = MM2ProjectPart1
if 1 == nargout
  internals = MM2ProjectInternals;
  return
end

% Specify the tolerances to be used in ODE45.
% OPTIONS should be passed as the final argument when calling ODE45.
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

% Vector of initial conditions for Task 1b
% Please arrange your eight initial values in the following order:
%  {xE(0), yE(0), uE(0), vE(0), <-- for Earth
%    xM(0), yM(0), uM(0), vM(0)} <-- for Mars
IC = [1.01671388;0;0;-6.176590574607843;1.666;0;0;-4.633];

% Run the code relating to Task 1b of Part 1
Task1_1b(options, IC);

% Run the code relating to Task 2 of Part 1
LaunchDates = 18:0.0025:19 %Array of launch dates for the period of 1 year.
Task1_2b(options, LaunchDates);

end

function [h, MinE, MinM, T, R] = Task1_1b(options, IC)
%TASK1_1B Computes nearest distances of Earth and Mars from the Sun.

if length(IC) ~= 8
  error('Array of initial conditions must have 8 components')
end

% Compute Earth-Mars orbits by calling EARTHMARSSYTEM function from ODE45
% Notes:
% (1) Run the simulation over an interval of four years.
% (2) Use the (array) IC for the initial conditions input argument.
% (3) Variable OPTIONS should be passed in as the final input argument.
timeSpan = [0,4]; %Setting start and end time in years.
[T,R] = ode45(@EarthMarsSystem, timeSpan, IC, options);

% Compute the closest distances of both Earth and Mars from the Sun
% during their orbits

%Calculating minimum distance of Earth from Sun
[noOfRRows, ~] = size(R); %Determining the number of rows of R
earthToSunDist = zeros(noOfRRows, 1); %Vector of Earth to Sun distances
for i = 1:noOfRRows
    earthToSunDist(i,1) = ((R(i,1))^2+(R(i,2))^2)^0.5;
end
MinE = min(earthToSunDist);

%Calculating minimum distance of Mars from Sun
marsToSunDist = zeros(noOfRRows, 1); %Vector of Mars to Sun distances
for i = 1:noOfRRows
    marsToSunDist(i,1) = ((R(i,5))^2+(R(i,6))^2)^0.5;
end
MinM = min(marsToSunDist);

% Plot Orbits
h = figure(1);
clf
hold on
plot(0, 0, 'oy', 'LineWidth', 7) % Sun's position
plot(R(:, 1), R(:, 2), 'b') % plot Earth's Orbit in blue
plot(R(:, 5), R(:, 6), 'r') % plot Mars' Orbit in red
xlabel('x')
ylabel('y')
title([
  'Orbits of Earth and Mars about the Sun: ' ...
  'Min. Earth Distance: ' num2str(MinE) '. ' ...
  'Min Mars Distance ' num2str(MinM)
  ])

end

function [h, OptLaunch, OptDistance, Reach] = Task1_2b(options, LaunchDates)
%TASK1_2 Computes maximum orbital distance and corresponding launch date.

% Generate figure in which to plot orbits
h = figure(2);
clf
hold on

% Pre-allocate storage for vector that contains orbital distances
Reach = zeros(1, length(LaunchDates));

% Iterate through the different launch dates
for j = 1 : length(LaunchDates)
  
  % Initialise the position of the Earth, Mars and probe by calling
  % provided function, RosettaIC.p, passing in each launch date in turn.
  % These twelve initial values are arranged in the following order:
  %  {xE(0), yE(0), uE(0), vE(0), <-- for Earth
  %    xM(0), yM(0), uM(0), vM(0), <-- for Mars
  %     xR(0), yR(0), uR(0), vR(0)} <-- for Rosetta2
  IC = Rosetta2IC(LaunchDates(j));
  
  % Determine the orbit of Earth, Mars and Rosetta2 by calling
  % EARTHMARSROSETTA2SYSTEM using ODE45.  
  % Notes:
  % (1) Run the simulation over an interval of five years.
  % (2) Use the (array) IC for the initial conditions input argument.
  % (3) Variable OPTIONS should be passed in as the final input argument.
  timeSpan = [0,5]; %Setting start and end time in years from launch date
  [T, R] = ode45(@EarthMarsRosetta2System, timeSpan, IC, options);
  
  %Determining the number of rows of R
  [noOfRRows, ~] = size(R);
  for i = 1:noOfRRows
      %Calculating the orbital distances of Rosetta from the Sun.
      rosettaOrbitDist(i) = ((R(i,9))^2+(R(i,10))^2)^0.5;
  end
  
  % Save max orbital distance of Rosetta2 for this launch date
  Reach(j) = max(rosettaOrbitDist);%Computing max. orbital distance of Rosetta2
  
  % Plot orbits for this launch date
  plot(0, 0, 'oy', 'LineWidth', 7) % Sun's position
  plot(R(:, 1), R(:, 2), '-b') % Plot Earth's orbit
  plot(R(:, 5), R(:, 6), '-r') % Plot Mars' orbit
  plot(R(:, 9), R(:, 10), '-g') % Plot Rosetta's orbit
  
end

% Compute the maximum orbital distance over these possible launch dates
% Note: You code may not look quite like this.
OptDistance = max(Reach); %Computing maximum orbital distance

%Calculating the optimal launch date.
%pos variable desribes the position in the Reach array that the max.
%orbital distance occurs. It also describes the position of the opt launch 
%date in the launch dates array.
[ ~ , pos] = max(Reach); %Tilde ignores the max. orbital dist. value.
OptLaunch = LaunchDates(pos); %Computing optimal launch date

% Annotate plot
title([
  'Optimal Probe Distance: ' num2str(OptDistance) '. ' ...
  'Optimal Launch Date: ' num2str(OptLaunch) ...
  ' years from 01/03/2012'
  ]);

end

function drdt = EarthMarsSystem(t, r)
%EARTHMARSSYSTEM Slopes for 1st order system corresponding to Equations 1.
% Derivatives should be stored in an order consistent with UNPACKSTATE.
% Note that the first argument is not actually required and may be replaced
% with a tilde ~.

%
% NB: Please use the following code if you get stuck:
% drdt = EarthMarsSystemReference(t, r) 
% return
%

%Loading the constant values
alphaConst = 1.272466364267360*10^-5;
betaConst = 39.65;
gammaConst = 1.1895*10^-4;

%Calculating the distance between earth and mars
rem = ( (r(5)-r(1))^2 + (r(6)-r(2))^2 )^0.5;
%Calculating the distance between earth and sun
rse = ((r(1))^2+(r(2))^2)^0.5;
%Calculating the distance between mars and sun
rsm = ((r(5))^2+(r(6))^2)^0.5;

%Setting up a column vector so that Matlab does not store values in a row
%vector.
drdt = zeros(8,1);

%Calculating the time derivatives
drdt(3) = -(alphaConst/(rem^3))*(r(1)-r(5)) - (betaConst*r(1))/(rse^3);
drdt(4) = -(alphaConst/(rem^3))*(r(2)-r(6)) - (betaConst*r(2))/(rse^3);
drdt(7) = -(gammaConst/(rem^3))*(r(5)-r(1)) - (betaConst*r(5))/(rsm^3);
drdt(8) = -(gammaConst/(rem^3))*(r(6)-r(2)) - (betaConst*r(6))/(rsm^3);

%Calculating the rest of the time derivatives
drdt(1) = r(3);
drdt(2) = r(4);
drdt(5) = r(7);
drdt(6) = r(8);
end

function drdt = EarthMarsRosetta2System(t, r)
%EARTHMARSROSETTA2SYSTEM Slopes for 1st order system corresponding to Equations 2.
% Derivatives are stored in an order consistent with UNPACKSTATE.

%Loading the constant values
alphaConst = 1.272466364267360*10^-5;
betaConst = 39.65;
gammaConst = 1.1895*10^-4;

%Calculating the distance between Rosetta2 probe and Mars
rrm = ((r(5)-r(9))^2+(r(6)-r(10))^2)^0.5;
%Calculating the distance between Rosetta2 probe and Earth
rre = ((r(1)-r(9))^2+(r(2)-r(10))^2)^0.5;
%Calculating the distance between Rosetta2 probe and the Sun
rrs = ((r(9))^2+(r(10))^2)^0.5;

%rosettaDerivative describes the slopes for only the motion of the 
%Rosetta2 probe.
rosettaDerivative = zeros(4,1); %Initialising variable.

%Calculating the time derivatives for the Rosetta2 probe
rosettaDerivative(1) = r(11);
rosettaDerivative(2) = r(12);
rosettaDerivative(3) = -(alphaConst/(rrm)^3)*(r(9)-r(5)) - (gammaConst/(rre)^3)*(r(9)-r(1)) - (betaConst*(r(9)))/((rrs)^3);
rosettaDerivative(4) = -(alphaConst/(rrm)^3)*(r(10)-r(6)) - (gammaConst/(rre)^3)*(r(10)-r(2)) - (betaConst*(r(10)))/((rrs)^3);

% HINT: Output order should consistent with the order of the initial conditions. This should be a vector containing 12 elements
drdt = [
  EarthMarsSystem(t, r); % a call to the function which encodes the EarthMars system, e.g. Equations 1, as built in Task 1.1b. This returns a vector containing 8 elements
  rosettaDerivative; % your new code to include the motion of the Rosetta probe, e.g. Equations 2.
  ];

end

function alpha = Alpha
%ALPHA The constant G*T^2*m_M/A_U^3
alpha = 1.272466364267360e-05;
end

function beta = Beta
%BETA The constant G*T^2*m_S/A_U^3.
beta = 39.65;
end

function gamma = Gamma
%GAMMA The constant G*T^2*m_E/A_U^3.
gamma = 1.1895e-04;
end

function normz = Norm(z)
%NORM Evalutae the norm/length of a vector.
%See also MATFUN/NORM.
normz = sqrt(z(1)^2 + z(2)^2);
end

function internals = MM2ProjectInternals
internals = struct( ...
  'EarthMarsSystem', @EarthMarsSystem, ...
  'EarthMarsRosetta2System', @EarthMarsRosetta2System, ...
  'Task1_1b', @Task1_1b, ...
  'Task1_2b', @Task1_2b);
end
