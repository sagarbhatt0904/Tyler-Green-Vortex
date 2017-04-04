function[x,y]=gridgen(N,st)
% Intializing the grid
if st==1
    dx=2*pi/(N-1);
    x1=0:dx:2*pi;
    y1=0:dx:2*pi;
else
    
    a1=(pi)*((st)^(-(N-1)/2));
    x1(1)=a1;
    y1(1)=a1;
    
    for w=2:((N-1)/2)
        
        x1(w)=a1*((st)^w);
        y1(w)=a1*((st)^w);
        
    end
    for w=((N+1)/2):N-1
        x1(w)=(2*pi)-x1(N-w);
        y1(w)=(2*pi)-x1(N-w);
    end
    x1(N)=2*pi;
    y1(N)=2*pi;
end
for i=1:N
    for j=1:N
        x(i,j)=x1(j);
        y(i,j)=y1(i);
    end
end
end
