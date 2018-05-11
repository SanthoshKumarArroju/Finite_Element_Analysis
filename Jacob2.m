function J2=Jacob2(nnel,dNdzi,dNdni,xcoord,ycoord);  % compute Jacobian (2-by-2 matrix)

J2 = zeros(2,2);

for i = 1: nnel
    
    J2(1,1) = J2(1,1) + (dNdzi(1,i)*xcoord(1,i));
    J2(1,2) = J2(1,2) + (dNdzi(1,i)*ycoord(1,i));
    J2(2,1) = J2(2,1) + (dNdni(1,i)*xcoord(1,i));
    J2(2,2) = J2(2,2) + (dNdni(1,i)*ycoord(1,i));
end 
end

