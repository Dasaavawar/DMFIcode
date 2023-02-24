escale=5000;

x1=cos(phi)*sint';
y1=sin(phi)*sint';
z1=ones(np,1)*cost';
x=zeros(np,nt,nr);
y=zeros(np,nt,nr);
z=zeros(np,nt,nr);
for i=1:nr
x(:,:,i)=r(i)*x1(:,:);
y(:,:,i)=r(i)*y1(:,:);
z(:,:,i)=r(i)*z1(:,:);
end


%%%%% SUBPLOT 1: TOP VIEW
sub1=subplot('position',[0 0 0.4 0.8])
daspect([1 1 1]);

innercore=surf(r(1)*x1(:,:),r(1)*y1(:,:),r(1)*z1(:,:),Br(:,:,1)/4)
set(innercore,'specularstrength',0.2);
hold on


% render the CMB radial field
cmbr=surf(r(nr)*x1(:,1:nt/2),r(nr)*y1(:,1:nt/2),r(nr)*z1(:,1:nt/2),Br(:,1:nt/2,nr))
caxis([-3,3])
alphad=(Br(:,1:nt/2,nr)).^2;
set(cmbr,'FaceAlpha','flat');
set(cmbr,'AlphaData',alphad);
set(cmbr,'specularstrength',0.2);
shading('flat')
alim([0 1])

if incr==1
brmin=min(min(Br(:,:,nr)));
step=0;
a=ones(np,nt);
while sum(sum(a))>30
step=step+0.01;
a=Br(:,:,nr)<(step*brmin);
end
mtheta=ones(np,1)*[1:nt];
mphi=[1:np]'*ones(1,nt);

sr=nr;
st=mtheta(a==1);
sp=mphi(a==1);
st=st(1:2:end);
sp=sp(1:2:end);

SX=r(sr)*sin(theta(st)).*cos(phi(sp));
SY=r(sr)*sin(theta(st)).*sin(phi(sp));
SZ=r(sr)*cos(theta(st));
end

verts=stream3(X,Y,Z,BX,BY,BZ,SX,SY,SZ);
verts2=stream3(X,Y,Z,-BX,-BY,-BZ,SX,SY,SZ);

draw=[verts verts2];
width=BX.^2+BY.^2+BZ.^2;

h=streamtube(draw,X,Y,Z,width/escale,[0 20]);

s=size(h);
s=s(1);
for i=1:s
XD=get(h(i),'XData');
YD=get(h(i),'YData');
ZD=get(h(i),'ZData');
BZD=interp3(X,Y,Z,BZ,XD,YD,ZD);
BZDmax=max(max(max(BZD)),-min(min(BZD)));
if BZDmax~=0 
ff=0.3*BZD./BZDmax;
else
ff=0;
end
ss=size(BZD);
CD=zeros(ss(1),ss(2),3);
CD(:,:,1)=0.7+ff;
CD(:,:,2)=0.7;
CD(:,:,3)=0.7-ff;
set(h(i),'CData',CD);
end

%set(h(:),'facecolor',[0.8,0.8,0.8])
set(h(:),'edgecolor','none')

for i=1:8
longitude=(i-1)*45/180*pi;
li=plot3(r(1)*sint*cos(longitude),r(1)*sint*sin(longitude),r(1)*cost)
set(li,'color','k')
     set(li,'linewidth',2)
end

for i=1:3
latitude=i*45/180*pi;
li=plot3(r(1)*sin(latitude)*cos(phi),r(1)*sin(latitude)*sin(phi),r(1)*cos(latitude)*ones(np,1))
set(li,'color','k')
set(li,'linewidth',2)
end


%%%% SUBPLOT 3: TOP VIEW INSERT
sub3=subplot('position',[0.4 0.8 0.1 0.2])
sc2=Br(:,:,nr);
spat_spec
Brsurf
spec_spat;
daspect([1 1 1]);

earthsurf=surf(r(1)*x1(:,:),r(1)*y1(:,:),r(1)*z1(:,:),150*sc2)
caxis([-3 3])
set(earthsurf,'specularstrength',0.2)
axis off
axis vis3d
shading flat
zoom(1.5)
view(0,90)
camlight(5,-10)
hold on

for i=1:8
longitude=(i-1)*45/180*pi;
li=plot3(r(1)*sint*cos(longitude),r(1)*sint*sin(longitude),r(1)*cost)
set(li,'color','k')
     set(li,'linewidth',2)
end

for i=1:3
latitude=i*45/180*pi;
li=plot3(r(1)*sin(latitude)*cos(phi),r(1)*sin(latitude)*sin(phi),r(1)*cos(latitude)*ones(np,1))
set(li,'color','k')
set(li,'linewidth',2)
end

%%%%%%% BACK TO COMPLETE SUBPLOT 1
set(sub1,'position',[0 0 0.5 1])
subplot(sub1)
zoom(1.5)
view(0,90)
camlight(5,-10)
axis off
axis vis3d

%%%%%%% SUBPLOT 2: SIDE VIEW
sub2=subplot('position',[0.50 0.0 0.39 0.79])
daspect([1 1 1]);

innercore=surf(r(1)*x1(:,:),r(1)*y1(:,:),r(1)*z1(:,:),Br(:,:,1)/4)
set(innercore,'specularstrength',0.2);
hold on


% render the CMB radial field
cmbr=surf(r(nr)*x1(np/2:np,:),r(nr)*y1(np/2:np,:),r(nr)*z1(np/2:np,:),Br(np/2:np,:,nr))
caxis([-3,3])
alphad=(Br(np/2:np,:,nr)).^2;
set(cmbr,'FaceAlpha','flat');
set(cmbr,'AlphaData',alphad);
set(cmbr,'specularstrength',0.2);
shading('flat')
alim([0 1])

h=streamtube(draw,X,Y,Z,width/escale,[0 20]);

s=size(h);
s=s(1);
for i=1:s
XD=get(h(i),'XData');
YD=get(h(i),'YData');
ZD=get(h(i),'ZData');
BZD=interp3(X,Y,Z,BZ,XD,YD,ZD);
BZDmax=max(max(max(BZD)),-min(min(BZD)));
if BZDmax~=0 
ff=0.3*BZD./BZDmax;
else
ff=0;
end
ss=size(BZD);
CD=zeros(ss(1),ss(2),3);
CD(:,:,1)=0.7+ff;
CD(:,:,2)=0.7;
CD(:,:,3)=0.7-ff;
set(h(i),'CData',CD);
end

set(h(:),'edgecolor','none')

for i=1:8
longitude=(i-1)*45/180*pi;
li=plot3(r(1)*sint*cos(longitude),r(1)*sint*sin(longitude),r(1)*cost)
set(li,'color','k')
     set(li,'linewidth',2)
end

for i=1:3
latitude=i*45/180*pi;
li=plot3(r(1)*sin(latitude)*cos(phi),r(1)*sin(latitude)*sin(phi),r(1)*cos(latitude)*ones(np,1))
set(li,'color','k')
set(li,'linewidth',2)
end


%%% DO THE DMFI

Exyz=BX.^2+BY.^2+BZ.^2;
for i=1:s/2
XD=get(h(i),'XData');
YD=get(h(i),'YData');
ZD=get(h(i),'ZData');
XD2=get(h(i+s/2),'XData');
YD2=get(h(i+s/2),'YData');
ZD2=get(h(i+s/2),'ZData');
ED=interp3(X,Y,Z,Exyz,XD,YD,ZD);
ED2=interp3(X,Y,Z,Exyz,XD2,YD2,ZD2);

[Exyzm,pos]=max(ED);
[Exyzm,pos2]=max(Exyzm);
[Exyzm2,pos3]=max(ED2);
[Exyzm2,pos4]=max(Exyzm2);
pos=pos(pos2);
pos3=pos3(pos4);

if Exyzm>Exyzm2
SX(i)=XD(pos,pos2);
SY(i)=YD(pos,pos2);
SZ(i)=ZD(pos,pos2);
else
SX(i)=XD2(pos3,pos4);
SY(i)=YD2(pos3,pos4);
SZ(i)=ZD2(pos3,pos4);
end

end

% SUBPLOT 4: SIDE VIEW INSERT

sub4=subplot('position',[0.9 0.8 0.1 0.2])
sc2=Br(:,:,nr);
spat_spec
Brsurf
spec_spat;
daspect([1 1 1]);

earthsurf=surf(r(1)*x1(:,:),r(1)*y1(:,:),r(1)*z1(:,:),150*sc2)
     caxis([-3 3])
     set(earthsurf,'specularstrength',0.2)
axis off
axis vis3d
shading flat
view(0,0)
zoom(1.5)
camlight(5,-10)
hold on

for i=1:8
longitude=(i-1)*45/180*pi;
li=plot3(r(1)*sint*cos(longitude),r(1)*sint*sin(longitude),r(1)*cost)
set(li,'color','k')
     set(li,'linewidth',2)
end

for i=1:3
latitude=i*45/180*pi;
li=plot3(r(1)*sin(latitude)*cos(phi),r(1)*sin(latitude)*sin(phi),r(1)*cos(latitude)*ones(np,1))
set(li,'color','k')
set(li,'linewidth',2)
end

% BACK TO COMPLETE SUBPLOT 2

set(sub2,'position',[0.50 0.0 0.49 0.99])
     subplot(sub2)
zoom(1.5)
view(0,0)
camlight(5,-10)
axis off
axis vis3d

