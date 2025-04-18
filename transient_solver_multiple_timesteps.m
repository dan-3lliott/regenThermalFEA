clc
clear all
close all

%import mesh data from Abaqus
nodes = importdata("nodes.txt");
nodes = nodes(:,2:end);
nodes = nodes .* 0.0254; %convert from inches to meters

elements = importdata("elements.txt");
elements = elements(:,2:end);

bartzNodes = importdata("bartzNodes.txt");
gnielinskiNodes = importdata("gnielinskiNodes.txt");

%define physical variables
thickness = 0.001; %m - will be used to go from 3D to 2D
Kxx = 325*thickness; %W/m*K
Kyy = 325*thickness; %W/m*K
hGnielinski = (2.13*10^4)*thickness; %W/m*K
hBartz = (1.00*10^4)*thickness; %W/m*K
rho = 8960*thickness; %kg/m^2
specificHeat = 385; %J/kg*K
coolantTemp = 444; %K
hotGasTemp = 2511; %K

%define transient variables
dTimes = [0.01 0.025 0.05 0.075 0.1]; %sec
timeDuration = 1; %sec
beta = 0.5; %Crank-Nicolson

%set up global matrices
nodeQty = size(nodes, 1);
elementQty = size(elements, 1);
K = zeros(nodeQty, nodeQty); %stiffness matrix
C = zeros(nodeQty, nodeQty); %specific heat matrix
D = zeros(nodeQty, 1) + 273.15; %temperature matrix, initial condition of ambient
avgTemp = [mean(D)];

%assemble global matrices
for i = 1:elementQty
    elementNodes = elements(i,:);
    %stiffness and specific heat
    [ke, ce] = elementLevelConstantMatrices(elementNodes, Kxx, Kyy, nodes, rho, specificHeat);
    K(elementNodes, elementNodes) = K(elementNodes, elementNodes) + ke;
    C(elementNodes, elementNodes) = C(elementNodes, elementNodes) + ce;
    %load
    R = globalLoadMatrix(nodeQty, nodes, gnielinskiNodes, hGnielinski, coolantTemp, bartzNodes, hBartz, hotGasTemp, D);
end

%start plotting average temperature over time
figure
%loop through all time steps
for dTime = dTimes
    %re-initialize R and D per initial conditions
    R = globalLoadMatrix(nodeQty, nodes, gnielinskiNodes, hGnielinski, coolantTemp, bartzNodes, hBartz, hotGasTemp, D);
    D = zeros(nodeQty, 1) + 273.15;
    avgTemp = [mean(D)];
    %re-define transient variables
    timeSteps = timeDuration/dTime;
    for i = 1:timeSteps
        Dprev = D;
        Rprev = R;
        %compute updated R based on temperatures and heat transfer coefficients
        R = globalLoadMatrix(nodeQty, nodes, gnielinskiNodes, hGnielinski, coolantTemp, bartzNodes, hBartz, hotGasTemp, Dprev);
        lhs = (C/dTime) + (beta * K);
        rhs = ((C/dTime) - ((1 - beta) * K)) * Dprev + (1 - beta)*Rprev + (beta*R);
        D = lhs\rhs;
        disp(i)
        avgTemp = [avgTemp mean(D)];
    end
    plot(0:dTime:timeDuration, avgTemp)
    hold on
end
grid on
xlabel('Time (sec)');
ylabel('Average Temperature (K)');
legend('dt = 0.01', 'dt = 0.025', 'dt = 0.05', 'dt = 0.075', 'dt = 0.1');
hold off

%function for element level matrices
function [k, c] = elementLevelConstantMatrices(elementNodes, Kxx, Kyy, nodes, rho, specificHeat)
    %stiffness
    elementCoords = nodes(elementNodes, :);
    xa = elementCoords(1,1);
    ya = elementCoords(1,2);
    xb = elementCoords(2,1);
    yb = elementCoords(2,2);
    xc = elementCoords(3,1);
    yc = elementCoords(3,2);
    A = 0.5*det([1 xa ya; 1 xb yb; 1 xc yc]);
    a1 = xb*yc - xc*yb;
    b1 = yb - yc;
    c1 = xc - xb;
    a2 = xc*ya - xa*yc;
    b2 = yc - ya;
    c2 = xa - xc;
    a3 = xa*yb - xb*ya;
    b3 = ya - yb;
    c3 = xb - xa;
    k = (1/(4*A)) .* [Kxx*(b1^2) + Kyy*(c1^2) Kxx*b1*b2 + Kyy*c1*c2   Kxx*b1*b3 + Kyy*c1*c3;
                     Kxx*b1*b2 + Kyy*c1*c2   Kxx*(b2^2) + Kyy*(c2^2) Kxx*b2*b3 + Kyy*c2*c3;
                     Kxx*b1*b3 + Kyy*c1*c3   Kxx*b2*b3 + Kyy*c2*c3   Kxx*(b3^2) + Kyy*(c3^2)];
    %specific heat 
    c = rho * specificHeat * (A/12) * [2 1 1; 1 2 1; 1 1 2];
end

%function for global load matrix
function [R] = globalLoadMatrix(nodeQty, nodes, gnielinskiNodes, hGnielinski, coolantTemp, bartzNodes, hBartz, hotGasTemp, Dprev) 
    R = zeros(nodeQty, 1);
    %gnielinski nodes
    for i = 1:(length(gnielinskiNodes)-1)
        n1 = gnielinskiNodes(i);
        n2 = gnielinskiNodes(i+1);
        le = norm(nodes(n1,:) - nodes(n2,:));
        T_avg = mean(Dprev([n1 n2]));
        %convert h to q
        qGniel = -hGnielinski * (T_avg - coolantTemp);
        R([n1 n2]) = R([n1 n2]) + (qGniel * le / 2) * [1; 1];
    end
    %bartz nodes
    for i = 1:(length(bartzNodes)-1)
        n1 = bartzNodes(i);
        n2 = bartzNodes(i+1);
        le = norm(nodes(n1,:) - nodes(n2,:));
        T_avg = mean(Dprev([n1 n2]));
        %convert h to q
        qBartz = -hBartz * (T_avg - hotGasTemp);
        R([n1 n2]) = R([n1 n2]) + (qBartz * le / 2) * [1; 1];
    end
end

%function for plotting the results
function plotMesh(nodes, elements, temperature)
    patch('Faces', elements, 'Vertices', nodes, 'FaceVertexCData', temperature, 'FaceColor', 'interp', 'EdgeColor', 'black');
    colormap(jet);
    a = colorbar;
    a.Label.String = 'Temperature (K)';
    clim([0 1000]);
    axis equal;
    xlabel('X (m)');
    ylabel('Y (m)');
    title('Regen Channel Thermal Analysis');
end