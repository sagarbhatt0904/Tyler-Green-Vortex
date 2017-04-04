function[u_new,v_new,p_new]=BC(N,u_new,v_new,p_new)
u_new(:,1)=0.5*(u_new(:,2)+u_new(:,N-1));
u_new(1,:)=0.5*(u_new(2,:)+u_new(N-1,:));
u_new(N,:)=u_new(1,:);
u_new(:,N)=u_new(:,1);

v_new(:,1)=0.5*(v_new(:,2)+v_new(:,N-1));
v_new(1,:)=0.5*(v_new(2,:)+v_new(N-1,:));
v_new(N,:)=v_new(1,:);
v_new(:,N)=v_new(:,1);

p_new(:,1)=0.5*(p_new(:,2)+p_new(:,N-1));
p_new(1,:)=0.5*(p_new(2,:)+p_new(N-1,:));
p_new(N,:)=p_new(1,:);
p_new(:,N)=p_new(:,1);

end
