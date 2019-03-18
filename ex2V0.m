clearvars
close all

h=4.5;
p=11e4;
F=3e5;
Y=2e11;
Area=250e-4;

nodes=0:h:4*h;
nodes=nodes';
elem=[1,2;2,3;3,4;4,5];
numNod=size(nodes, 1);
numElem=size(elem,1);

Ke=Y*Area/h*[1,-1;-1,1];
K = zeros(numNod);
Q = zeros(numNod,1);
u = zeros(numNod,1);


fixedNodes=1;
freeNodes=setdiff(1:numNod,fixedNodes);

Ke=Area*Y/h*[1,-1;-1,1];

for e=1:numElem
    row=[elem(e,1),elem(e,2)];
    col=row;
    if e==2
        K(row,col)=K(row,col)+(5e9*125e-4)/h*[1,-1;-1,1];
    else
        K(row, col)=K(row, col) + Ke;
    end
end

%Natural B.C
Q(freeNodes)=-2*p;
Q(numNod)=-F;

%Essential B.C
u(fixedNodes)=0;


Qm=Q(freeNodes)-K(freeNodes,fixedNodes)*u(fixedNodes);
Km=K(freeNodes, freeNodes);

um=Km\Qm;
u(freeNodes)=um;
format short e
u


%PostProcess
rectionForces=K*u-Q

for e=1:numElem
    displ(e)=u(elem(e,2))-u(elem(e,1));
    stress(e)=Y*displ(e)/h;
    force(e)=Area*stress(e);
end

display(['displ', 'force', 'stress'])
[displ', force', stress']