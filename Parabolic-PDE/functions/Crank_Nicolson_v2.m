function u_new=Crank_Nicolson_v2(u_old,t_alpha)
    n=length(u_old);
    N=n-2;
    
    alpha=t_alpha*ones(1,N)/2;
    beta=(1+t_alpha)*ones(1,N);
    gamma=alpha;

    alpha(1)=0;
    gamma(N)=0;

    v=zeros(1,N);
    w=zeros(1,N);
    delta=zeros(1,N);
    

   
        for i=1:N
            delta(i)=alpha(i)*u_old(i)+(1-t_alpha)*u_old(i+1)+gamma(i)*u_old(i+2);
        end

    v(1)=gamma(1)/beta(1);
    w(1)=delta(1)/beta(1);


        for i=2:N
            v(i)=gamma(i)/(beta(i)-alpha(i)*v(i-1));
            w(i)=(delta(i)+alpha(i)*w(i-1))/(beta(i)-alpha(i)*v(i-1));
        end

    u_new(N)=w(N);

        for i=N-1:-1:1
            u_new(i)=v(i)*u_new(i+1)+w(i);
        end


end
