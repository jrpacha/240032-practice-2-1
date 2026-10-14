clearvars
close all

H = 3.0;           %m
h = 1.0;           %m
P = 2.0e3;         %N
E = 28.0e9;        %N / m^2
w = 25e3;          %N / m^3
R = 24.0e3 * 0.25; %N Força a l'extrem superior de la columna 

%nodes = (0:h:H)';        %posició dels nodes

nodes = linspace(0,H,4)'; %alternativament amb la funció linspace
elem = [1, 2;             %matriu de connectivitat 
    2, 3;
    3, 4];

alphaK = -E / 6; betaK = 3*E / 4; 
alphaF = w / 6; betaF  = -3*w / 4; 

numNods = size(nodes,1);
numElem = size(elem, 1);

K = zeros(numNods);
F = zeros(numNods, 1);
Q = zeros(numNods, 1);
u = zeros(numNods, 1);

%Ke = A* E * [1, -1; 
%    -1, 1] / h;

%Fe = [0; 0];

for e = 1:numElem
    rows = [elem(e,1), elem(e,2)];
    cols = rows;

    x1 = nodes(elem(e,1)); x2 = nodes(elem(e,2));
    
    %Càlcul de Ke
    Ke = (0.5 * (x1 + x2) * alphaK + betaK) * [1, -1; -1, 1] / h;
 
    %Càlcul de Fe
    Fe = (alphaF*[2*x1 + x2; x1 + 2*x2] / 6 + betaF * [1; 1] /2) * h;

    %Acoplem les matrius Ke i els vectors Fe de cada element
    K(rows, cols) = K(rows, cols) + Ke
    F(rows) = F(rows) + Fe;     
end

%Boundary Conditions (B.C.)
fixedNods = 1;
freeNods = setdiff(1:numNods,fixedNods);

%Natural B.C.
Q(3) = -2.0*P;
Q(4) = -R;

%Essential B.C.
u(1) = 0.0;

%Reduced system
Km = K(freeNods,freeNods);
Fm = F(freeNods) + Q(freeNods) - K(freeNods,fixedNods)*u(fixedNods);

%Solve the reduced system
um = Km \ Fm;
u(freeNods) = um;

%Post Process: reaction forces
Q = K * u - F;

format short e
fprintf("%11s%12s%11s\n","nodePos (m)","U (mm)","Q (N)")
fprintf("%7.3f%17.4e%11.3e\n",[nodes,u*1000, Q]')