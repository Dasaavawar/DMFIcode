%    Routine to load magic data into matlab
%    Works for binary data !
%    ok so let's scan the data first
%    Julien Aubert 11/02, jaubert@gwdg.de
    
     fid=fopen(fname);     
     version=char(fread(fid,32,'char')');
     tit=char(fread(fid,64,'char')');
     gridpar=fread(fid,9,'float');
     nr=gridpar(4)
     nt=gridpar(5)
     np=gridpar(6)
     nr_ic=gridpar(7);
     azsym=gridpar(8)
     nblock=gridpar(9);


     Vr=zeros(np,nt,nr);
     Vp=zeros(np,nt,nr);
     Vt=zeros(np,nt,nr);
     Br=zeros(np,nt,nr+nr_ic);
     Bp=zeros(np,nt,nr+nr_ic);
     Bt=zeros(np,nt,nr+nr_ic);
     T=zeros(np,nt,nr);

     phypar=fread(fid,6,'float');
     Ek=phypar(2)
     Pr=phypar(3)
     Ra=phypar(1)
     Pm=phypar(4)

     dummy=fread(fid,2,'float');
     theta=fread(fid,nt,'float');
     cost=cos(theta);
     sint=sin(theta);

     r=zeros(nr+nr_ic,1);

     for ir=1:nr*nblock
     dummy=fread(fid,2,'float');
     raddata=fread(fid,4,'float');
     rlevel=raddata(1);
     r(rlevel+1)=raddata(2);
     
     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,np,'float');
     T(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,np,'float');
     Vr(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,np,'float');
     Vt(:,i,rlevel+1)=dummy;
     end
     
     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,np,'float');
     Vp(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,np,'float');
     Br(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,np,'float');
     Bt(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,np,'float');
     Bp(:,i,rlevel+1)=dummy;
     end

     end

     for ir=1:nr_ic

     dummy=fread(fid,2,'float');
     raddata=fread(fid,4,'float');
     rlevel=raddata(1);
     r(rlevel+1)=raddata(2);

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,np,'float');
     Br(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,np,'float');
     Bt(:,i,rlevel+1)=dummy;
     end

     for i=raddata(3):raddata(4)
     dummy=fread(fid,2,'float');
     dummy=fread(fid,np,'float');
     Bp(:,i,rlevel+1)=dummy;
     end

     end

     fclose(fid);

%    now for sorting the data out of this strange hemispherical
%    way that the magic code uses

     for ir=1:nr
     dummy2=T(:,:,ir);
     for i=1:nt/2
     T(:,i,ir)=dummy2(:,2*(i-1)+1);
     T(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end

     dummy2=Vr(:,:,ir);
     for i=1:nt/2
     Vr(:,i,ir)=dummy2(:,2*(i-1)+1);
     Vr(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end

     dummy2=Vp(:,:,ir);
     for i=1:nt/2
     Vp(:,i,ir)=dummy2(:,2*(i-1)+1);
     Vp(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end
     
     dummy2=Vt(:,:,ir);
     for i=1:nt/2
     Vt(:,i,ir)=dummy2(:,2*(i-1)+1);
     Vt(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end
     end

     for ir=1:nr+nr_ic

     dummy2=Br(:,:,ir);
     for i=1:nt/2
     Br(:,i,ir)=dummy2(:,2*(i-1)+1);
     Br(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end
     
     dummy2=Bt(:,:,ir);
     for i=1:nt/2
     Bt(:,i,ir)=dummy2(:,2*(i-1)+1);
     Bt(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end

     dummy2=Bp(:,:,ir);
     for i=1:nt/2
     Bp(:,i,ir)=dummy2(:,2*(i-1)+1);
     Bp(:,i+nt/2,ir)=dummy2(:,nt-1-2*(i-1)+1);
     end

     end

% RENORMALIZE r!

   r=r/(1-0.35);     

% Phi

   phi=[1:np]'/np*2*pi/azsym;

