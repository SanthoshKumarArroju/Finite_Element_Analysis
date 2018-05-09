
function [Ns,Ns_elements,dNdzi,dNdni]=dNs_Q4(xi,ni); % compute shape functions and derivatives at sampling point

Ns = (1/4)*[(1-xi)*(1-ni) (1+xi)*(1-ni) (1+xi)*(1+ni) (1-xi)*(1+ni)];

Ns_elements = [Ns(1),0,Ns(3),0;
               0,Ns(2), 0, Ns(4);
               Ns(3),0,Ns(1),0;
               0,Ns(4), 0, Ns(2)];

dNdzi =(1/4)*[-(1-ni) (1-ni) (1+ni) -(1+ni)];

dNdni =(1/4)*[-(1-xi) -(1+xi) (1+xi) (1-xi)];

end