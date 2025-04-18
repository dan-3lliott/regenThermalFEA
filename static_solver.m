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
Kxx = 325; %W/mK
Kyy = 325; %W/mK
qGnielinski = -(1.8*10^7); %W/m^2
tHotWall = 850; %K

%set up global matrices
nodeQty = size(nodes, 1);
elementQty = size(elements, 1);
K = zeros(nodeQty, nodeQty);
R = zeros(nodeQty, 1);

%assemble global matrices
for i = 1:elementQty
    elementNodes = elements(i,:);
    [ke, re] = elementLevelMatrices(elementNodes, Kxx, Kyy, nodes, gnielinskiNodes, qGnielinski);
    K(elementNodes, elementNodes) = K(elementNodes, elementNodes) + ke;
    R(elementNodes) = R(elementNodes) + re;
end

%enforce hot wall temperature boundary condition
for i = 1:length(bartzNodes)
   K(bartzNodes(i),:) = 0;
   K(bartzNodes(i),bartzNodes(i)) = 1;
   R(bartzNodes(i)) = tHotWall;
end

%compute temperature matrix
D = K\R;

%plot results
plotMesh(nodes, elements, D)

%function for element level matrices
function [k, r] = elementLevelMatrices(elementNodes, Kxx, Kyy, nodes, gnielinskiNodes, qGnielinski)
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
    %load
    r = zeros(3,1);
    [gnielMask, ~] = ismember(elementNodes, gnielinskiNodes);
    if (sum(gnielMask) > 1)
        xe = elementCoords(gnielMask,1);
        ye = elementCoords(gnielMask,2);
        le = sqrt((xe(1)-xe(2))^2 + (ye(1)-ye(2))^2);
        r(gnielMask) = r(gnielMask) + (qGnielinski * le/2);
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