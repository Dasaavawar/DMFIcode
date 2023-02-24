mcost=(ones(nr+nr_ic,1)*cost')';
msint=(ones(nr+nr_ic,1)*sint')';

mmcost=zeros(np,nt,nr+nr_ic);
mmsint=zeros(np,nt,nr+nr_ic);
on=ones(nt,nr+nr_ic)  ;     
       

for i=1:np
mmcost(i,:,:)=mcost;
mmsint(i,:,:)=msint;
mmcosphi(i,:,:)=on*cos(phi(i));
mmsinphi(i,:,:)=on*sin(phi(i));
end

a=[-1:0.02:1]*r(1);
b=[-1:0.02:1]*r(1);
c=[-1:0.02:1]*r(1);

[X Y Z]=meshgrid(a,b,c);

[PHI LAT R]=cart2sph(X,Y,Z);
TH=pi/2-LAT;
PHI(PHI<0)=2*pi+PHI(PHI<0);

for i=1:101
i
for j=1:101
for k=1:101
[a,phim(i,j,k)]=min(abs(PHI(i,j,k)-phi));
[b,thm(i,j,k)]=min(abs(TH(i,j,k)-theta));
[c,rm(i,j,k)]=min(abs(R(i,j,k)-r));

cond(i,j,k)=sqrt(X(i,j,k).^2+Y(i,j,k).^2+Z(i,j,k).^2)<r(1)*1.01;

end
end
end
