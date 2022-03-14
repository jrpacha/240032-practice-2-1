clearvars
close all

%Consider the BVP
%
% -((x+3)'u)' + u  = x, 0 < x < L
%            u'(0) = du0,
%             u(L) = uL
%
% by the FEM, taking N=100 linear elements, compute the nodal solutions 
% u(k), k=1,2,...,N+1 taking L=50, du0=4.94 and uL= 1.0, and give its 
% maximum value, max{u(k),k=1,2,...,N+1}.
%
% Hint: The value of the FEM approximation for u at node 51 is
% 2.4895614e+01

%Data
numDiv=100;    %number of divisions
L=50.0;        %0 <= x <= L.   
alpha=1.0;
beta=3.0;      %a1(x)=alpha*x+beta
gamma=1.0;    %a0(x)=gamma
a=1.0;         %f(x)=a*x;

du0 = 4.94;    %value of u' at x = 0. Natural BC u'(0) = 4.94, 
uL = 1.0;      %Essential B.C. u(L)=1.0;

%Exercise:
%(a) Show that for a1(x)=alpha*x+beta, and linar elements (2 nodes), the
%K1e term of the stiffness matrix of the elemenet e is given by 
%     Ke1 = alpha*(x1e+x2e)*[1,-1;-1,1]/(2*he) + beta*[1,-1;-1,1]/he,
%where x1e and x2e are the positions of the 1st and the 2nd node of the
%element, respectively, and he is its length.
%
%(b) Show that, for f(x)=a*x+b, and linear elements (2 nodes) the load 
%vector of element e is given by
%    Fe = a*he*[2*x1e+x2e;x1e+2*x2e]/6 + b*he*[1;1]/2,
%where x1e and x2e are the positions of the 1st and the 2nd node of the
%element, respectively, and he is its length.

%-------------------------------------------------------------------------
%Geometry: nodes & elements
X0=0.0; X1=L;
h=(X1-X0)/numDiv;
nodes=(X0:h:X1)';
elem=[1:numDiv;2:numDiv+1]';
numNodes=size(nodes,1);
numElem=size(elem,1);
%-------------------------------------------------------------------------

%Assembly
Q=zeros(numNodes,1);
F=zeros(numNodes,1);
u=zeros(numNodes,1);
K=zeros(numNodes);

Ke0=beta*[1,-1;-1,1]/h + gamma*h*[2 1;1 2]/6;  %This part of the stiffness
                                               %matrix is the same for all
                                               %the elements
for e=1:numElem
    rows=[elem(e,1),elem(e,2)]; cols=rows;
    x1=nodes(rows(1,1)); x2=nodes(rows(1,2));
    Ke=Ke0+0.5*alpha*(x1+x2)*[1,-1;-1,1]/h;
    Fe=a*h*[2*x1+x2;x1+2*x2]/6;
    K(rows,cols)=K(rows,cols)+Ke;
    F(rows)=F(rows)+Fe;
end

%BC
fixedNods=numNodes;
freeNods=setdiff(1:numNodes,fixedNods);

%Natural B.C.
Q(1) = -(alpha*nodes(1)+beta)*du0; %The other are Q(2)=...=Q(numNod-1)=0

%Essential B.C
u(fixedNods)=uL;

%Set the reduced system
Qm = Q(freeNods) + F(freeNods) - K(freeNods,fixedNods)*u(fixedNods);
Km = K(freeNods,freeNods);
um = Km\Qm;
u(freeNods)=um;

%Solution
fprintf('Solution:\n')
fprintf('max{u(k),k=1....,N+1} = %.4e\n',max(u));
fprintf('Hint. value of the FEM approx. for\n u at node 51, u(51) = %.7e\n',...
    u(51))