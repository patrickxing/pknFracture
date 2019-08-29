%% code by Pengju, to solve the PKN model without fluid leakoff
% Solve numerically by finite difference method
%------------------------------------------------------------
%% Define the parameters
clc
clear
E=100e6;
nu=0.48;
Ep=E/(1-nu^2);
mu=0.3; % fluid viscosity
G=E/(2*(1+nu));
Q=1e-7;
q1=Q/2;
h=0.05;
B=pi*G/256/(1-nu)/mu;
t=50.001;
L=0.68*(G*q1^3/(1-nu)/mu/h^4)^(0.2)*t.^0.8; % In term of G
W0t=2.5*((1-nu)*mu*q1^2/G/h)^0.2*t.^0.2; % In term of G
% W0tp=2.5*(2*mu*q1^2/Ep/h)^0.2*t.^0.2;  % In term of Ep
% Lp=0.68*(Ep*q1^3/2/mu/h^4)^(0.2)*t^0.8; % in term of Ep
%% Theoretical solution of W with time and x
% Get the fitting relation from Nordgren's paper
Wxt_W0t=[1	0.9677762	0.9265587	0.88760835	0.8418864	0.7871481	0.72792804	0.65291995	0.56442124	0.41955724	0.3065856		0];
x_L=[0	0.100821756	0.20108797	0.2991283	0.39938703	0.5007457	0.59875244	0.70007753	0.79803586	0.9048194	0.9581441		1.00];
N=100;
x_L_p=linspace(0,1,N+1);  %End point of the element
ii=1;
for ii=1:N
    x_L_pM(ii)=(x_L_p(ii)+x_L_p(ii+1))/2; % Dimensionless middle point of the element
end
Wxt_W0tIp=interp1(x_L,Wxt_W0t,x_L_pM); %interpolation results of the crack width for the middle point
Wxt_PK=4*((1-nu)*mu*q1*L/pi/G*(1-x_L_pM)).^0.25; %PK solution
% figure (1)
% plot(x_L_pM,Wxt_W0tIp) % The curve obtained 
% hold on
% plot(x_L,Wxt_W0t,'r') % Dimensionless curve from Nordgren's paper, obtained by Plot digitizer

% %The comparison between dimensional form 
% %Comparison between numerical solution and PK solution
figure (2)
Wxt=W0t*Wxt_W0tIp;
plot(x_L_pM*L,Wxt) %Nordgren Numerical solution
hold on 
plot(x_L_pM*L,Wxt_PK,'r')  %PK solution
%% Start the numerical method
% Initial condition, m=0, N=1
W(1)=0;
dx=0.01; % assume the element length is 0.01
dT0=0.001; % assume the initial time increase is 0.010s, note that the initial value should be small
dW(1)=4*dT0*q1/(pi*dx*h);
W=[W+dW 0];
% check the volume balance
dV=dx*pi*h*dW/4; % calculated volume based on updated crack width
errorR(1)=(dV-q1*dT0)/dV; % the error is zero for the intial step

% Start the calculation, the outer is time loop
M=100; % M is the number of the time loops, the number of elements is N=M+1
dT=ones(M,1)*0.5; % For PKN model, the crack length is proportional to t^(4/5), almost linearly
for m=1:M
    N=m+1; % N is the number of elements
    % Initialize the matrix
    A=zeros(N,N);
    C=zeros(N,1);
    %while (abs(error)>0.05) % this loop is calculating the crack increase
    %=========================================================================
    %Principly, we should use while loop to check the global volume balance for each time loop
    %However, due to the amolst linear relation of time vs length in PKN model, the global balance
    %error is very small provided we choose proper element length and time step. 
    %So there is no need for the while loop
    %=========================================================================
        % The loop over the elements
        for n=1:N
            if n==1 % the first row of the matirx
                A(1,1)=-4*W(1)^3-pi*h/4*2*dx^2/dT(m)/B;
                A(1,2)=4*W(2)^3;
                C(1)=-2*W(2)^4+2*W(1)^4-2*q1*dx/B;
            elseif n==N % last row of the matrix
                A(N,N)=-pi*h/4*2*dx^2/dT(m)/B;
                A(N,N-1)=4*W(N-1)^3;
                %C(N)=3*W(N)^4-2*W(N-1)^4;
                C(N)=-2*W(N-1)^4;
            end
            if N>2 && n~=1 && n~=N
                A(n,n+1)=4*W(n+1)^3;
                A(n,n)=-pi*h/4*2*dx^2/dT(m)/B-8*W(n)^3;
                A(n,n-1)=4*W(n-1)^3;
                C(n)=-2*W(n+1)^4+4*W(n)^4-2*W(n-1)^4;
            end
            dW=(A\C)'; 
        end
        if min(dW)<0 % the variation of crack width should be positive
            break
        end
        dVC=sum(dx*pi*h*dW/4); % calculated volume change based on updated crack width
        dVP=q1*dT(m); % volume change based on pump rate
        error=(dVC-dVP)/dVC; % volume variation error 
        errorR(m+1)=error; % array of the volume variation error
        W=[W+dW 0]; % At the end of m step, the number of elements is m+2
        VC=sum(W)*dx*pi*h/4; % calculated total volume based on updated crack width
        VP=q1*(sum(dT(1:m))+dT0); % total volume based on pump rate
        errorS=(VC-VP)/VC; % total volume error 
        errorSum(m)=errorS;% array total volume error 
    %end
end
% Make sure the volume balance error not be large
errorRMax=max(abs(errorR));
errorSumMax=max(abs(errorSum));

% output the numerical crack width change with x at the end of simulation
% compare the numerical solution with the solution of Fig.4 in Nordgren's paper (dimensional form)
xNum=0.01:0.01:(M+2)*0.01;
hold on
plot(xNum,W, 'ro')   
legend('Nordgren solution', 'Numerical solution of this code')
xlabel('x (m)')
ylabel('crack width (m)')

