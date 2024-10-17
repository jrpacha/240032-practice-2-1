%B 
% 
% 
% 
% 
% 
% asic 
% Quiz P2-1
%Compute the value at the point x=1 for the solution of the diff. equation
%1.2 u''+ 0.2u =0, for 0<=x<=2 with boundary conditions u(0)=1.2, u'(2)=0
%Use 10 linear elements to approximate the solution.

%Trieu-ne una:

%(a) 1.60806680e+00 (*)
%(b) 1.14303688e+00
%(c) 1.61514427e+00
%(d) 3.21397712e+00

clearvars;
close all

a1 = -1.2;      %Ull!!!! Signe!!! l'equació model és 
                % -(a1(x) u')' + a0(x) u = f(x)
a0 = 0.2;

xini = 0;
xfin = 2;
numDiv = 10; 

uIni = 1.2;     %B.C. at x = xini (essential B.C.!!!)
dudxFin = 0.0;  %B.C. at x = xfin (QFin = a1(xfin) * dudxFin,
                %                  natural B.C.!!!)

xp = 1.0;       %Point where the solution is interpolated at

%Compute the solution
numElem = numDiv;
h = (xfin - xini)/numElem;

nodes = (xini:h:xfin)'; %node coordinates (as column vector)

elem = zeros(numElem,2);

for i=1:numElem 
       elem(i,:) = [i, i+1];
end	

numNod = size(nodes,1);

numElem = size(elem,1);

Ke = a1/h*[1,-1;-1,1] + a0*h/6*[2,1;1,2]; %local stiff matrix (constant!!)

K = zeros(numNod); %initialize the global Stiff Matrix
F = zeros(numNod,1); %initialize the internal forces global vector
Q = zeros(numNod,1); %initialize the global Q vector

for i=1:numElem 
       rows = [elem(i,1); elem(i,2)];
       colums = rows;	
       K(rows,colums) = K(rows,colums)+Ke; %assembly Stiff matrix
       % Assembly of the  non-constant F terms
       %nod1 = elem(i,1);
       %nod2 = elem(i,2);
       %Fe = 0;
       %F(rows,1) = F(rows,1)+Fe;
end
Fini = F; %copy the original F values for the postprocess

%fixedNodes = [1,numNod]; %Fixed nodes, here the first and the last (global num.)
%                         %Nooooooo!!!

%Ull!!!! Només està fixada la u en l'extrem de l'esquerra (x = 0). L'altra 
%condició de contorn ens fixa la derivada en l'extrem de la dreta (x = 2)
%Per tant, en aquest cas les condicions de contorn són u(1) = 1.2 (u(0) = 
%1.2) en l'equació i Q(end) = a1*0 (Q(xfin) = a1(xfin)*u'(xfin) en l'equació). 
% Recorda que quan ens donen la derivades podem calcular les Q's!

fixedNodes = 1; %La u nomé està fixada al 1er node!!!!

freeNodes = setdiff(1:numNod,fixedNodes); %Complementary of fixed nodes
% modify the linear system, here BC are NOT 0.

%Natural B.C. (fix Q's)
Q(end) = a1*dudxFin;  % **En aquest cas** no caldria, ja que el 
                      % vector duDxFin = 0, i el vector Q està 
                      %inicialitzat a zero... 

%Essential B.C. (fix u's)
u = zeros(numNod,1);  %initialize u vector
u(1) = uIni;          %BC at x=xini

%u(numNod) = 0; %BC at x=xfin
Q = Q-K(:,fixedNodes)*u(fixedNodes);  
 
Km = K(freeNodes,freeNodes);

Fm = F(freeNodes)+Q(freeNodes);

% solve the System
format long e; %just to a better view of the numbers

um = Km\Fm;
u(freeNodes) = um;

%interp1(u,1) % Ull!!!!! a la funció interp1 li has de passar també la
              % posició dels nodes...
interpU = interp1(nodes, u, 1);
%interpU

%Fancy outptut: don't do this at exams.
fprintf("Solution: interpolated value of u at x = %.2f: %.8e\n",...
    xp,interpU)