%  And now for the field tubes!!!!
%  First compute cartesian field components
%    Julien Aubert 11/02, jaubert@gwdg.de

%  Interpolate B on the mesh. Ok, nearest neighbor is crude
%  but it does the job....

Bs=Br.*mmsint+Bt.*mmcost;
Bx=Bs.*mmcosphi-Bp.*mmsinphi;
By=Bs.*mmsinphi+Bp.*mmcosphi;
Bz=Br.*mmcost-Bt.*mmsint;

for i=1:101
for j=1:101
for k=1:101


BX(i,j,k)=Bx(phim(i,j,k),thm(i,j,k),rm(i,j,k))*cond(i,j,k);
BY(i,j,k)=By(phim(i,j,k),thm(i,j,k),rm(i,j,k))*cond(i,j,k);
BZ(i,j,k)=Bz(phim(i,j,k),thm(i,j,k),rm(i,j,k))*cond(i,j,k);


end
end
end
