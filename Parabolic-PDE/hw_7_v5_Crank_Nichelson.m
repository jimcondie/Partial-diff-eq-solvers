%PARABOLIC PDE SOLVER 
%using implict Crank_Nicolson scheme
%version 5- plots error analysis
%v4 uses 'improved' Crank_Nicolson_v2 function in an attempt to
%improve accuracy

clc
clear

%set approximation parameters
t_len=4000;
x_len=51;  %keep odd to included midpoint


%set spatial boundaries
lower_x_b=0;
upper_x_b=1;

%set time bounds
lower_t_b=0;
upper_t_b=.5;



t=linspace(lower_t_b,upper_t_b,t_len);
x=linspace(lower_x_b,upper_x_b,x_len);

t_step=t(2)-t(1); %dimensionless
x_step=x(2)-x(1); %dimensionless
alpha=t_step/x_step^2;  %stable for all values

u=zeros(x_len,t_len);
analytic=zeros(x_len,t_len);
temp_error=zeros(1,x_len);
error=zeros(1,t_len);
maxiter=100; %for series expansion limit in analytic soltn

%set initial conditions here
for i=1:x_len
    if x(i)<.5
         u(i,1)=2*x(i);
    else
         u(i,1)=2*(1-x(i));
    end
end


%estimate forward in time 
for j=1:(t_len-1)
    
    %set BC's here
    u(1,j)=0;
    u(x_len,j)=0;
    
    %set the rest of u at for one advancement in time
    u(2:x_len-1,j+1)=Crank_Nicolson_v2(u(:,j),alpha);

        %find analytic solution
        for index=1:x_len
            for a=1:maxiter
                analytic(index,j)=analytic(index,j)+1/a^2*sin(a*pi/2)*sin(a*pi*x(index))*exp(-a^2*pi^2*t(j));
            end
        end
        
     analytic(:,j)=8/pi^2*analytic(:,j);
     
     %calculate error and store maximum error in x as function of t
     temp_error=abs((analytic(:,j)'-u(:,j)')/u(:,j)')*100; %percent error wrt analytical soltn
     error(j)=max(temp_error);
     
    %plot results for certain times
        if ( j==1 || mod(j,500)==0) 
                           
        
        figure
        plot(x,u(:,j),x,analytic(:,j),'--')
         ylim([0,1])
             if(j==1)
             legend({'Crank-Nicolson approx.','analytic soltn.'},'Location','South')
             else
             legend({'Crank-Nicolson approx.','analytic soltn.'},'Location','Northeast')
             end
        ylabel('dimensionless variable u')
        xlabel('dimensionless variable in space')
        title(['U at t = ',num2str(t(j)) ,' dimensionless variable in time'])
        end
end


%error analysis
figure
plot(t,error)
% ylim([3e-6,4e-4])
title(['Crank Nicolson error with t-step: ',num2str(t_step),' and x-step: ',num2str(x_step)])
xlabel('t, dimensionless variable in time')
ylabel('error % (wrt analytical solution)')
