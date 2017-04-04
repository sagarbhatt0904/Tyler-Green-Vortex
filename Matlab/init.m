function[xvel,yvel,Press,xvel1,yvel1,Press1]=init(N,Re,t,t1,x,y)
%Initializing U,V, Pressure
for i=1:N
            for j=1:N
                xvel(i,j)=-exp((-2*t/Re))*cos(x(1,i))*sin(y(j,1));
                yvel(i,j)=exp((-2*t)/Re)*cos(y(j,1))*sin(x(1,i));
                Press(i,j)=-((cos(2*x(1,i))+cos(2*y(j,1)))*exp((-4*t)/Re))/4;
                xvel1(i,j)=-exp((-2*t1)/Re)*cos(x(1,i))*sin(y(j,1));
                yvel1(i,j)=exp((-2*t1)/Re)*cos(y(j,1))*sin(x(1,i));
                Press1(i,j)=-((cos(2*x(1,i))+cos(2*y(j,1)))*exp((-4*t1)/Re))/4;
                
            end
        end
end
