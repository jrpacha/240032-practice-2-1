clearvars
close all

%Compute the value at the point x=1 for the solution of the diff. equation 
%1.7 u''+ 0.4u =0, for 0<=x<=2 with boundary conditions u(0)=1.8, u'(2)=0.
%Use 10 linear elements to approximate the solution.

%Data
numDiv=10;      %number of divisions
L=2.0           %  0 <= x <= 2.   
a1 = -1.7;
a0 = 0.4;
u0 = 1.8;       %value of u at x = 0. Essential BC u(0) = 1.8, 
dudx2 = 0.0;    %Natural B.C. u'(2) = 0;
%-------------------------------------------------------------------------
%Geometry: nodes & elements
h=L/numDiv;
nodes=0:h:L;
nodes=nodes(:);
%-------------------------------------------------------------------------

elem=zeros(numDiv,2);
for i=1:numDiv
    elem(i,:)=[i,i+1];
end

numNod=size(nodes,1);
numElem=size(elem,1);

%Assembly
Q=zeros(numNod,1);
u=zeros(numNod,1);
K=zeros(numNod);

Ke = a1*[1, -1;-1,1]/h + a0*h*[2 1; 1 2]/6.0; %This Ke is the same for 
                                              %for all the elements, since
                                              %coefficients are constants.
                                              %See T3-MN-FEM2D.pdf,page 56

for e=1:numElem
    %Ke=localStiffnessMatrix1D(E,A,nodes,elem,e);
    rows=[elem(e,1),elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke;
end

u = zeros(numNod,1);

%Natural B.C.
Q(numNod) = a1*dudx2; 

%Essential B.C
fixedNodes=1;
freeNodes=setdiff(1:numNod,fixedNodes);
u(fixedNodes)=u0;

%Set the reduced system
Qm = Q(freeNodes) - K(freeNodes,fixedNodes)*u(fixedNodes);
Km = K(freeNodes,freeNodes);
um = Km\Qm;
u(freeNodes)=um;
format long e
%u
[nodes,u]
u(6)