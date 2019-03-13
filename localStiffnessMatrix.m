function Ke = localStiffnessMatrix(E,A,nodes,elem,e)
L = abs(nodes(elem(e,2))-nodes(elem(e,1)));
Ke=A(e)*E(e)/L*[1,-1;-1,1];
end
