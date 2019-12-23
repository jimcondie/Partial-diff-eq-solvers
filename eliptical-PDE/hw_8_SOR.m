%hw 8
%solving eliptical PDE using Successive Over Relaxation (SOR)

clc 
clear

%numerical parameters
steps=400;
max_iter=100000;
max_error=1e-12;
figures=10;
omega=1.97;  %stable for 1<omega<2

%set dimensions
size=3; %cm;

%set initial guess
guess=1;  %initial constant guess
phi_new=guess*ones(steps,steps);

%set initial conditions

%top border at 5 volts

phi_new(steps,:)=5*ones(1,steps);

%other borders are at zero 
phi_new(1,:)=zeros(1,steps);
phi_new(2:steps-1,1)=zeros(steps-2,1);
phi_new(2:steps-1,steps)=zeros(steps-2,1);

%create copy to keep track of old iterations
phi_old=phi_new;


%Gauss-Seidel guesses for new iterations
done=false;
for iter=1:max_iter
    if(done==true)
        break
    end
    
    for i=2:steps-1
        for j=2:steps-1
            phi_new(j,i)=(1-omega)*phi_old(j,i)+0.25*omega*(phi_new(j,i-1)+phi_old(j,i+1)+phi_new(j-1,i)+phi_old(j+1,i));
        end
    end
    
    %check for convergence
    error_matrix=abs(phi_new-phi_old)./phi_old;
    error=max(max(error_matrix));
    if (error<max_error)
        done=true;
    end
    
    phi_old=phi_new;

    %plot various progressions
    
    %plots just final result
    if(done==true)
        
    %plots less often for higher order iterations
     %if(iter==1 || (iter<10 && mod(iter,2)==0)||(iter<100 && mod(iter,10)==0)||( iter<1000 && mod(iter,100)==0) || mod(iter, max_iter/figures)==0 || done==true)
        x=linspace(0,size,steps);
        y=linspace(0,size,steps);
     [X,Y]=meshgrid(x,y);
    
     figure
    C=phi_new;
    surf(X,Y,phi_new,C,'FaceAlpha',.79,'EdgeColor', 'none','FaceColor','interp')
    view(-48,25)
    c=colorbar;
    c.Label.String = 'Volts';
    camlight left
    axis tight
        if(done==true)
        title(['final SOR iteration: ',num2str(iter),', mesh size:',num2str(steps),'^2 steps', ', initial guess phi=: ',num2str(guess)]) 
        else
        title(['SOR iteration: ',num2str(iter)])
        end
    zlabel({'Volts'});
    ylabel({'cm'});
    xlabel({'cm'});
    end
end


