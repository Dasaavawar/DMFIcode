%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                       %%
%%      PPPPP         A        RRRRR       OOOOO    DDDD    Y     Y      %%
%%      P    P       A A       R    R     O     O   D   D    Y   Y       %%
%%      P    P      A   A      R   R      O     O   D    D    Y Y        %%
%%      PPPP       A     A     RRRR       O     O   D    D     Y         %%
%%      P          AAAAAAA     R  R       O     O   D    D     Y         %%
%%      P         A       A    R   R      O     O   D   D      Y         %%
%%      P         A       A    R    R      OOOOO    DDDD       Y         %%
%%                                                                       %%
%%                                                                       %%
%%                        PArallel  ROtating  DYnamos                    %%
%%                        --        --        --                         %%
%%                                                                       %%
%%                                                                       %%
%%  Julien Aubert, Emmanuel Dormy                                        %%
%%  Philippe Cardin                                                      %%
%%  Patrick Stoclet                                                      %%
%%  Vincent Morin                                                        %%
%%  Laure Goudard                                                        %%
%%                                                                       %%
%%  Version JA-1.0 of May 2006                                           %%
%%                                                                       %%
%%  This is a matlab routine to load parody data                         %%
%%  Type fname='your_G_file_name'                                        %%
%%  and  parodyload                                                      %%
%%                                                                       %%
%%                                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(fname);
  
dummy=fread(fid,1,'float');
version=fread(fid,1,'float');

dummy=fread(fid,2,'float');
phypar=fread(fid,10,'float')';
time=phypar(1);
DeltaU=phypar(2);
Coriolis=phypar(3);
Lorentz=phypar(4);
Buoyancy=phypar(5);
ForcingU=phypar(6);
DeltaT=phypar(7);
ForcingT=phypar(8);
DeltaB=phypar(9);
ForcingB=phypar(10);

Ek=1/Coriolis;
Ra=Buoyancy*Ek;
Pm=1/DeltaB;
Pr=1/DeltaT;

dummy=fread(fid,2,'float');
gridpar=fread(fid,4,'float');
nr=gridpar(1);
nt=gridpar(2);
np=gridpar(3);
azsym=gridpar(4);

dummy=fread(fid,2,'float');
r=fread(fid,nr,'float');

dummy=fread(fid,2,'float');
theta=fread(fid,nt,'float');
sint=sin(theta);
cost=cos(theta);
phi=[1:np]'/np*2*pi/azsym;

dummy=fread(fid,2,'float');

Vr=zeros(np,nt,nr);
Vt=zeros(np,nt,nr);
Vp=zeros(np,nt,nr);
Br=zeros(np,nt,nr);
Bt=zeros(np,nt,nr);
Bp=zeros(np,nt,nr);
T=zeros(np,nt,nr);

for ir=1:nr
for it=1:nt
Vr(:,it,ir)=fread(fid,np,'float');
dummy=fread(fid,2,'float');
Vt(:,it,ir)=fread(fid,np,'float');
dummy=fread(fid,2,'float');
Vp(:,it,ir)=fread(fid,np,'float');
dummy=fread(fid,2,'float');
Br(:,it,ir)=fread(fid,np,'float');
dummy=fread(fid,2,'float');
Bt(:,it,ir)=fread(fid,np,'float');
dummy=fread(fid,2,'float');
Bp(:,it,ir)=fread(fid,np,'float');
dummy=fread(fid,2,'float');
T (:,it,ir)=fread(fid,np,'float');
dummy=fread(fid,2,'float');
end
end

fclose(fid);
