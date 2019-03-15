function Ke = localStiffnessMatrix1D(E,A,nodes,elem,e)
%LOCALSTIFFNESSMATRIX1D
%Compute the local Stiffness matrix for the element e in 1D
%FEM.
%INPUT
%     E: row array holding the Young moduli of the elements.
%     A: row array holding the section area of each element.  
% nodes: column array with the coordinates if the nodes.  
%  elem: connectivity matrix. It has 2 columns and as many rows
%        as elements.
%     e: number of the current element for which the local
%        stiffness matrix is computed.
%OUTPUT
%    Ke: local stiffness matrix for the element e
%    
L = abs(nodes(elem(e,2))-nodes(elem(e,1)));
Ke=A(e)*E(e)/L*[1,-1;-1,1];
end
