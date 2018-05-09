clc; clear all; close all;
% Q4 Element formation scalar field problem  
 
%% To read node locations in local coordinate system
Filename_node= 'NLIST.lis';
Fid_node = fopen(Filename_node,'r');
input_flag_string_node = '   NODE        X           Y           Z         THXY    THYZ    THZX';

while ~feof(Fid_node)
    
   Node_input = fgetl(Fid_node);
   flag_input = strfind(Node_input, input_flag_string_node)
   
   if flag_input
       
       for i = 1:515       % give last node number in NLIST
           
           SS = fgetl(Fid_node);
           SS = str2num(SS);
           NN(i,:) = SS;

       end
   end
    
end
    fclose(Fid_node)
    X = NN(:,2);
    Y = NN(:,3);
    Z = NN(:,4);
    THETA_Z = NN(:,5);
    THETA_X = NN(:,6);
    THETA_Y = NN(:,7);
    
    gcoord = [X,Y]

%% To read elements from the input file    
Filename_element= 'ELIST.lis';
Fid_element = fopen(Filename_element,'r');

input_flag_string_element = '    ELEM MAT TYP REL ESY SEC        NODES';

while ~feof(Fid_element)
    
   element_input = fgetl(Fid_element);
   flag_input = strfind(element_input, input_flag_string_element)
   
   if flag_input
       
       for i = 1:408       % give last node number in NLIST
           
           SS = fgetl(Fid_element);          
           SS = str2num(SS);
           EE(i,:) = SS;

       end
   end
    
end
fclose(Fid_element)

nodes = EE(:,7:end)

nodes_dof = nodes;

%% plot the model
   gcoord_x = gcoord(1:end,1);
   gcoord_y = gcoord(1:end,2);  
   X_node = gcoord_x(nodes(1:end,1:end))';
   Y_node = gcoord_y(nodes(1:end,1:end))';
%    figure
   plot(X_node,Y_node)
   axis([-5,105,-20,20])
    
%%  FE model parameter information

nel=length(nodes(:,1))                 % number of elements
nnel=length(nodes(1,:))                 % number of nodes per element
ndof=1                                  % number of dofs per node           scalar field problem
nnode=length(gcoord(:,1))               % total number of nodes in system
ngp=3                                  % use 2x2 integration rule
sdof=nnode*ndof                        % total system dofs  
edof=nnel*ndof                          % dofs per element

%%  Physical properties  
rho = 2700        % density of material
A = 0.9655       % elemental area   units: mm^2
Ee = 2.1e5     % young modulus units: N/mm^2
v = 0.33       % poissons ratio

Dm = Ee/(1-v^2)*[1 v ; v 1];          % for plane stress


Fg=zeros(sdof,1);             % initialization of system force vector
Kg=zeros(sdof,sdof);          % initialization of system matrix
Mg=zeros(sdof,sdof);
elm_dofs=zeros(nnel*ndof,1);  % initialization of index vector
[gp,gpw]=GP_quadrature1D(ngp);  % sampling points & weights


for iel=1:nel             %   1:nel loop for the total number of elements

    for i=1:nnel
        nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
        xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
        ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
    end
    
    Ke=zeros(edof,edof);         % initialization of element Stiffness matrix to zero
    Mass_e = zeros(edof,edof);   % initialization of element Mass matrix to zero
    
  
    for int_xi=1:ngp
        xi=gp(int_xi);                  % sampling point in x-axis
        wtxi=gpw(int_xi);               % weight in x-axis
        for int_ni=1:ngp
            ni=gp(int_ni);              % sampling point in y-axis
            wtni=gpw(int_ni) ;          % weight in y-axis

            [Ns,Ns_elements,dNdzi,dNdni]=dNs_Q4(xi,ni); % compute shape functions and derivatives at sampling point
            J2=Jacob2(nnel,dNdzi,dNdni,xcoord,ycoord);  % compute Jacobian
            detJ2=det(J2);                 % determinant of Jacobian
            invJ2=inv(J2);                 % inverse of Jacobian matrix
            Bm=invJ2*[dNdzi;dNdni]         % compute the strain_displacement B matrix
            
            Ke=Ke+Bm'*Dm*Bm*detJ2*1*wtxi*wtni % unit thickness
            Mass_e = Mass_e + rho*A*Ns'*Ns*1*detJ2*wtxi*wtni;
        
        end
       
    end  % end of numerical integration loop
    
     elm_dofs = nodes_dof(iel,:);
     
     Kg(elm_dofs,elm_dofs)=Kg(elm_dofs,elm_dofs)+Ke;  
     
     Mg(elm_dofs,elm_dofs)=Mg(elm_dofs,elm_dofs)+Mass_e;
             
     
end   % end of element loops


Kg   % global stiffness matrix
Mg   % global mass matrix

alpha = 10; beta = 2/alpha;
Cg = alpha*Mg + beta*Kg   % global Damping matrix, using proportional damping method

%% Loads and Boundary Conditions
Fg(2)=150;Fg(102:107)=150;

Kg_beforeBC=Kg;
Mg_beforeBC=Mg;
Cg_beforeBC=Cg;

bcdof=[130 131 132 133];         % dof only one direction scalar field 
bcval=[0 0 0 0];                % Boundary conditions

for i=1:length(bcdof)
    bci=bcdof(i);
    bcv=bcval(i);
    Kg(bci,:)=Kg(bci,:)*0; Kg(bci,bci)=1, Fg(bci)=bcv
    Mg(bci,:)=Mg(bci,:)*0; Mg(bci,bci)=1
    Cg(bci,:)=Cg(bci,:)*0; Cg(bci,bci)=1
end

Kg_afterBC = Kg;
Mg_afterBC = Mg;
Cg_afterBC = Cg;

%% Solution displacement
Phi=inv(Kg_afterBC)*Fg;

Phi_x = Phi(1:2:end,:);  % odd row matrix, x-direction solution
Phi_y = Phi(2:2:end,:);  % even row matrix, y-direction solution
 
  X_nodalsol = Phi_x(nodes(1,1:end))';
  Y_nodalsol = Phi_y(nodes(1,1:end))';
  
%   figure
%   plot(X_nodalsol,Y_nodalsol)            % * denote node points


%% Stress

% for iel=1:nel             %   1:nel loop for the total number of elements

%     for i=1:nnel
%         nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
%         xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
%         ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
%     end
%  
%     elm_dofs = nodes_dof(iel,:)
%     
%     for int_xi=1:ngp
%         xi=gp(int_xi);                  % sampling point in x-axis
%         wtxi=gpw(int_xi);               % weight in x-axis
%         for int_ni=1:ngp
%             ni=gp(int_ni);              % sampling point in y-axis
%             wtni=gpw(int_ni) ;          % weight in y-axis
% 
%             [Ns,Ns_elements,dNdzi,dNdni,dNdzini]=dNs_Q8(xi,ni); % compute shape functions and derivatives at sampling point
% 
%             Bm=[1 0 0 0; 0 0 0 1; 0 1 1 0]*[invJ2 zeros(2); zeros(2) invJ2]*[dNdzi;dNdni;dNdzj;dNdnj]         % compute the strain_displacement B matrix
%           
% %             edof = 1:nel;
% %                for iel = 1:nel
% %          elm_dofs = nodes_dof(iel,:)
%          Phi(elm_dofs,1)
%          strain_e = Bm*Phi(elm_dofs,1)
%          Phi(nodes_dof(2,:),1)
%          for ii = 1:24;
%          strain_total = zeros(3,nel)
%          strain_total(3,iel) = strain_total(3,iel) + strain_e


%         end      
%     end  % end of numerical integration loop

%     strain_total = zeros(3*nel,1);
    
%      elm_dofs = nodess(iel,:);
%      
% %     elm_dofs=Elem_dofs(nd,nnel,ndof);  % extract element dofs for the Kg assembly
%      Kg(elm_dofs,elm_dofs)=Kg(elm_dofs,elm_dofs)+Ke;  
%      
%      elm_dofs = nodess(iel,:);
% 
%      Mg(elm_dofs,elm_dofs)=Mg(elm_dofs,elm_dofs)+Mass_e;
     
    
% end   % end of element loops


%% Potential energy
Strain_energy = 1/2 * Phi'*Kg*Phi
workdone = Phi'*Fg
Potential_energy = Strain_energy - workdone 


%% Eigen Values and Eigen Vectors
  
  L_forMg = chol(Mg_beforeBC);
  L_forKg = chol(Kg_beforeBC);
  L_forCg = chol(Cg_beforeBC);

    Mgt = inv(L_forMg)*Mg_beforeBC*inv(L_forMg');
    Kgt = inv(L_forKg)*Kg_beforeBC*inv(L_forKg');
    Cgt = inv(L_forCg)*Cg_beforeBC*inv(L_forCg');
    
   [eigVector,eigValues]=eig(Kgt,Mgt)

   Wn = sqrt(diag(eigValues));     % natural frequency rad/sec, eigen values
   mode_shape = real(eigVector(:,1));   % mode shape,  eigen vectors
   
%    mmode = mode_shape(nodes(1:10,1:end))
   
   

%    t = 0:0.1:10; wd = 50;  f = Fg;    %Fg = f*cos(wd*t);        
%    
% X = (Kg_beforeBC - Mg_beforeBC*wd^2)\f;
%% using vibration tool box

% [v_vector,w_natrualfreq,zeta_dampingratio]=vtb4_3(Mgt, Cgt ,Kgt)   


%% Un-damped Homogeneous Vibration     at Fundamental Frequency Wn
Xo = 10; Vo = 25;
t = 0:0.5:514;
X_undamp = (sqrt(Xo^2*abs(Wn(1)).^2 + Vo^2)/abs(Wn(1))) * sin(t*abs(Wn(1)) + atan( (Xo*abs(Wn(1))) / Vo ))
figure
plot(t,X_undamp)
title('Un Damped Homogeneous, Displacement vs time plot')
xlabel('Time - sec'), ylabel('Displacement - mm'), grid on


%% Damped Homogeneous Vibration
Xo = 10; Vo = 25;                 % initial conditions
zeta = Cgt / 2*sqrt(Mgt*Kgt);
Wd = sqrt(1-zeta^2)*Wn;

A = sqrt( Xo^2 +(( Vo + zeta(1)*abs(Wn(1))*Xo) /abs(Wd(1)) )^2 );
Phi_angle = atan(Xo*abs(Wd(1)) + (Vo + zeta*abs(Wn(1))*Xo));

t = 0:0.01:514;
X_A = A* exp(-zeta_dampingratio(1)*abs(Wn(1))*t) ;
X_damped = X_A .* (sin(abs(Wd(1))*t + Phi_angle(1)));
figure
plot(t(1:250),X_damped(1:250))
title('Damped Homogeneous Vibrations, Displacement vs time plot')
xlabel('Time - sec'), ylabel('Displacement - mm'), grid on

%% Impulse forced vibration
dt = 0.01; F = 150;
F_cap = F*dt
t = 0:0.01:514;

AA = (F_cap /(Mgt(1)*abs(Wd(1)) )* exp(-zeta_dampingratio(1)*abs(Wn(1))*t))
X_impulse = AA .*sin(abs(Wd(1))*t)
figure
plot(t(1:300),X_impulse(1:300),'b')
title('Impulse forced vibration, Displacement vs time plot')
xlabel('Time - sec'), ylabel('Displacement - mm'), grid on

%% FFT Frequency Domain Analysis
figure
H_undamped = fft(X_undamp);
Freq = 0:0.5:514;
plot(Freq,20*log10(abs(H_undamped)))
title('undamped vibration in Frequency Domain, FRF vs Frequency plot')
xlabel('Frequency - Hz'), ylabel('FRF'), grid on


figure
H_damped = fft(X_damped);
Freq = 0:0.01:514;
plot(Freq,20*log10(abs(H_damped)))
title('Damped vibration in Frequency Domain, FRF vs Frequency plot')
xlabel('Frequency - Hz'), ylabel('FRF'), grid on


figure
H_impulse = fft(X_impulse);
Freq = 0:0.01:514;
plot(Freq,20*log10(abs(H_impulse)))
title('Impulse forced Damped vibration in Frequency Domain, FRF vs Frequency plot')
xlabel('Frequency - Hz'), ylabel('FRF'), grid on

