clearvars
close all

h=4.5;      %m
P=11.0e4;   %N
F=3.0e5;    %N
Y=2.0e11;   %N/m^2
Area=250.0e-4; %m^2

nodes=linspace(0,4*h,5);
nodes=nodes';
elem=[1,2;
      2,3;
      3,4;
      4,5];

numNod=size(nodes,1);
numElem=size(elem,1);

E=Y*ones(1,numElem);
A=Area*ones(1,numElem);

u=zeros(numNod,1);
Q=zeros(numNod,1);

K=zeros(numNod);
for e=1:numElem
    Ke=localStiffnessMatrix(E,A,nodes,elem,e);
    rows=[elem(e,1),elem(e,2)];
    files=rows;
    K(files,rows)=K(files,rows)+Ke;
end

fixedNod=1;
freeNod=setdiff(1:numNod,fixedNod);

%Natural B.C.
Q(freeNod)=-2*P;
Q(numNod)=-F;

%Essential B.C
u(fixedNod)=0.0;

Qm = Q(freeNod) - K(freeNod,fixedNod)*u(fixedNod);
Km = K(freeNod,freeNod);
um = Km\Qm;
u(freeNod)=um;
format short e
u

%Post-process
reactForces = K*u-Q

displ=zeros(numElem,1);
stress=zeros(numElem,1);
force=zeros(numElem,1);

for e=1:numElem
    displ(e) = u(elem(e,2))-u(elem(e,1));
    L = abs(nodes(elem(e,2))-nodes(elem(e,1)));
    stress(e) = E(e)*displ(e)/L;
    force(e) = A(e)*stress(e);
end

fprintf("%12s%14s%14s\n",'displ.','force','stress')
fprintf('%14.4e%14.4e%14.4e\n',[displ,force,stress]')



