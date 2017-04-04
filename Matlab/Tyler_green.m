clear all;
close all;
N1=[41];                          %Grid size
st1=[1]; 			  %stretching ratio
Re=100;
eps=0.01;                         %Epsilon for dissipation
L=30;
t=0;
t1=0.1;
Etd=0;
Es=0;
Es_st=0;
dt=[0.005];			  %Real time step	
for o=1:1
    for m=1:1
        for l=1:1
        N=N1(m);
        st=st1(o);
        [x,y]=gridgen(N,st);  % Grid generation
        [xvel,yvel,Press,xvel1,yvel1,Press1]=init(N,Re,t,t1,x,y); % Initial conditions
        [JC,zx,ey,zy,ex]=metric(N,x,y);   %Metric Calculation
        u_new1=xvel;
        u_new2=xvel;
        u_new3=xvel;
        u_new=xvel;
        u_k1=xvel;
        u_k2=xvel;
        u_k3=xvel;
        u_k=xvel;
        u_n1=xvel;
        u_n2=xvel;
        u_n3=xvel;
        u_n=xvel;
        u_old1=xvel;
        u_old2=xvel;
        u_old3=xvel;
        u_old=xvel;
        v_new1=yvel;
        v_new2=yvel;
        v_new3=yvel;
        v_new=yvel;
        v_k1=yvel;
        v_k2=yvel;
        v_k3=yvel;
        v_k=yvel;
        v_n1=yvel;
        v_n2=yvel;
        v_n3=yvel;
        v_n=yvel;
        v_old1=yvel;
        v_old2=yvel;
        v_old3=yvel;
        v_old=yvel;
        p_new1=Press;
        p_new2=Press;
        p_new3=Press;
        p_new=Press;
        p_n=Press;
        dj=1;
        
            for ti=0:dt(l):t1
                ps=0;
                while dj>=10^(-5)
                    for i=1:N
                        for j=1:N
                            dtau(i,j)=0.00005;
                        end
                    end
                    
                    [rcs,rus,rvs]=RHS(N,JC,u_k,v_k,Press,zx,ey,ex,zy,Re,ep,eps);
                    RHSu=(((-3*u_k+4*u_n-u_old)/(2*dt(l)))+rus);
                    RHSv=(((-3*v_k+4*v_n-v_old)/(2*dt(l)))+rvs);
                    for i=1:N
                        for j=1:N
                            p_new1=Press+0.25*(dtau(i,j)*rcs(i,j));
                            u_new1=u_k+0.25*(dtau(i,j)*RHSu(i,j));
                            v_new1=v_k+0.25*(dtau(i,j)*RHSv(i,j));
                        end
                    end
                    [u_new1,v_new1,p_new1]=BC(N,u_new1,v_new1,p_new1);
                    
                    [RHSc1,RHSu1,RHSv1]=RHS(N,JC,u_new1,v_new1,p_new1,zx,ey,ex,zy,Re,ep,eps);
                    RHSu11=(((-3*u_k1+4*u_n1-u_old1)/(2*dt(l)))+RHSu1);
                    RHSv11=(((-3*v_k1+4*v_n1-v_old1)/(2*dt(l)))+RHSv1);
                    for i=1:N
                        for j=1:N
                            p_new2(i,j)=Press(i,j)+0.33*(dtau(i,j)*RHSc1(i,j));
                            u_new2(i,j)=u_k(i,j)+0.33*(dtau(i,j)*RHSu11(i,j));
                            v_new2(i,j)=v_k(i,j)+0.33*(dtau(i,j)*RHSv11(i,j));
                        end
                    end
                    [u_new2,v_new2,p_new2]=BC(N,u_new2,v_new2,p_new2);
                    
                    [RHSc2,RHSu2,RHSv2]=RHS(N,JC,u_new2,v_new2,p_new2,zx,ey,ex,zy,Re,ep,eps);
                    RHSu22=(((-3*u_k2+4*u_n2-u_old2)/(2*dt(l)))+RHSu2);
                    RHSv22=(((-3*v_k2+4*v_n2-v_old2)/(2*dt(l)))+RHSv2);
                    for i=1:N
                        for j=1:N
                            p_new3(i,j)=Press(i,j)+0.5*(dtau(i,j)*RHSc2(i,j));
                            u_new3(i,j)=u_k(i,j)+0.5*(dtau(i,j)*RHSu22(i,j));
                            v_new3(i,j)=v_k(i,j)+0.5*(dtau(i,j)*RHSv22(i,j));
                        end
                    end
                    [u_new3,v_new3,p_new3]=BC(N,u_new3,v_new3,p_new3);
                    
                    [RHSc3,RHSu3,RHSv3]=RHS(N,JC,u_new3,v_new3,p_new3,zx,ey,ex,zy,Re,ep,eps);
                    RHSu33=(((-3*u_k3+4*u_n3-u_old3)/(2*dt(l)))+RHSu3);
                    RHSv33=(((-3*v_k3+4*v_n3-v_old3)/(2*dt(l)))+RHSv3);
                    for i=1:N
                        for j=1:N
                            p_new(i,j)=Press(i,j)+(dtau(i,j)*RHSc3(i,j));
                            u_new(i,j)=u_k(i,j)+(dtau(i,j)*RHSu33(i,j));
                            v_new(i,j)=v_k(i,j)+(dtau(i,j)*RHSv33(i,j));
                        end
                    end
                    [u_new,v_new,p_new]=BC(N,u_new,v_new,p_new);
                    
                    
                    dj= (norm((xvel1-u_new)));
                    u_old1=u_n1;
                    u_n1=u_k1;
                    u_k1=u_new1;
                    u_old2=u_n2;
                    u_n2=u_k2;
                    u_k2=u_new2;
                    u_old3=u_n3;
                    u_n3=u_k3;
                    u_k3=u_new3;
                    u_old=u_n;
                    u_n=u_k;
                    u_k=u_new;
                    v_old1=v_n1;
                    v_n1=v_k1;
                    v_k1=v_new1;
                    v_old2=v_n2;
                    v_n2=v_k2;
                    v_k2=v_new2;
                    v_old3=v_n3;
                    v_n3=v_k3;
                    v_k3=v_new3;
                    v_old=v_n;
                    v_n=v_k;
                    v_k=v_new;
                    
                    p1=p_new1;
                    p2=p_new2;
                    p3=p_new3;
                    Press=p_new;                    
                end  
            end
%	  Computing temporal error
                        for i=1:N
                            for j=1:N
                                Et(i,j)=sqrt(((Press1(i,j)-p_new(i,j))^2)+((xvel1(i,j)-u_new(i,j))^2)+((yvel1(i,j)-v_new(i,j))^2));
                                Etd=Et(i,j)+Etd;
                            end
                        end          
            
                        Err_t(l)=((N)^(-2))*Etd;
                        Etd=0;
        end
%         Computing spatial error
                 for i=1:N
                     for j=1:N
                         E(i,j)=(((Press1(i,j)-p_new(i,j))^2)+((xvel1(i,j)-u_new(i,j))^2)+((yvel1(i,j)-v_new(i,j))^2));
                         Es=E(i,j)+Es;
                     end
                 end
                 Err(m)=((N)^(-2))*Es;
                 Es=0;
    end
%             Computing spatial error for stretching
                for i=1:N
                    for j=1:N
                        E(i,j)=(((Press1(i,j)-p_new(i,j))^2)+((xvel1(i,j)-u_new(i,j))^2)+((yvel1(i,j)-v_new(i,j))^2));
                        Es_st=E(i,j)+Es_st;
                    end
                end
                Err_st(o)=((N)^(-2))*Es_st;
                Es_st=0;

end

