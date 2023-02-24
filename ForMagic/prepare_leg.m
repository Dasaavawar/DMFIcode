% This routines prepares the gauss points and fully normalized legendre polynomials
% Julien Aubert Summer 2005 adpated from routines in MAGIC


% First get the truncation
nalias=20;
n_phi_max=np;
n_theta_max=nt;
minc=azsym;
l_max=fix((nalias*n_theta_max)/30);
	  m_max=fix((l_max/minc))*minc;
lm_max=fix(m_max*(l_max+1)/minc-m_max*(m_max-minc)/(2*minc)+(l_max+1-m_max));

% compute gauss points and weights

dpi=double(pi);
N=double(n_theta_max);
M=double((N+1)/2);
XXM=0.d0;
XXL=1.d0;
EPS=3.D-16;

for I=1:M
ZZ=double(cos( dpi*( (double(I)-0.25D0)/(double(N)+0.5D0)) ));

ZZ1=0;
while (double(abs(ZZ-ZZ1)>EPS))

         P1=1.D0;
         P2=0.D0;
         for J=1:N
	 P3=P2;
	 P2=P1;
	 P1=( double(2*J-1)*ZZ*P2-double(J-1)*P3 )/double(J);
	 end

         PP=double(N)*(ZZ*P1-P2)/(ZZ*ZZ-1.D0);
         ZZ1=ZZ;
         ZZ=ZZ1-P1/PP;
end

         XX(I)=double(acos(XXM+XXL*ZZ));
XX(N+1-I)=double(acos(XXM-XXL*ZZ));
W(I)=2.D0*XXL/((1.D0-ZZ*ZZ)*PP*PP);
W(N+1-I)=W(I);

end

theta=XX';

% compute the legendre polynomials

m0=minc;
max_degree=l_max+1;
max_order=m_max;

for it=1:nt
colat=XX(it);

pos=0;
for m=0:m0:max_order;
if m==0
dnorm=1;
else
dnorm=sqrt(2);
end

fac=1.d0;
for j=3:2:2*m+1;
fac=fac*double(j)/double(j-1);
end;

plm=sqrt(fac);
if (sin(colat)~=0.d0) 
plm=plm*(-sin(colat))^m;
elseif( m~=0 ) then
plm=0.d0;
end
 
l=m;
pos=pos+1;
plma(it,pos) = dnorm*plm;
plm1=0.d0;
 

for l=m+1:max_degree;
plm2=plm1;
plm1=plm;
plm= cos(colat)* sqrt( double( (2*l-1)*(2*l+1) ) / double( (l-m)*(l+m) )  ) .* plm1 - sqrt( double( (2*l+1)*(l+m-1)*(l-m-1) ) / double( (2*l-3)*(l-m)*(l+m) ) ) .* plm2;
pos=pos+1;
plma(it,pos) = dnorm*plm;
end;

% additional plm(max_degree+1) stored in dtheta_plma(it,pos) to save space...
l=max_degree+1;
plm2=plm1;
plm1=plm;
plm= cos(colat)* sqrt( double( (2*l-1)*(2*l+1) ) / double( (l-m)*(l+m) )  ) * plm1 - sqrt( double( (2*l+1)*(l+m-1)*(l-m-1) ) / double( (2*l-3)*(l-m)*(l+m) ) ) * plm2;
dtheta_plma(it,pos)=dnorm*plm;
end;

% evaluation of sin(theta)*theta deriv of plm using a recurrence


% l=m contribution
pos=0;
for m=0:m0:max_order
if m==0
dnorm=1;
else
dnorm=sqrt(2);
end

l=m;
pos=pos+1;
dtheta_plma(it,pos)= l/sqrt(double(2*l+3)) * plma(it,pos+1);

% l>m contribution
for l=m+1:max_degree-1
pos=pos+1;
dtheta_plma(it,pos)= l*sqrt( double((l+m+1)*(l-m+1)) / double((2*l+1)*(2*l+3)) ) * plma(it,pos+1)  - (l+1)*sqrt( double((l+m)*(l-m)) / double((2*l-1)*(2*l+1))) * plma(it,pos-1);
end;

% l=max_degree contribution, note usage of dtheta_plma(it,pos) instead of plm(max_degree+1)
l=max_degree;
pos=pos+1;
dtheta_plma(it,pos)= l*sqrt( double((l+m+1)*(l-m+1)) / double((2*l+1)*(2*l+3)) ) * dtheta_plma(it,pos)  - (l+1)*sqrt( double((l+m)*(l-m)) / double((2*l-1)*(2*l+1)) ) * plma(it,pos-1);

% end loop over order
end;

% end loop over colatitudes
end;




for n_theta=1:n_theta_max/2
lm=0;
lmp=0;
for m=0:minc:m_max;
for l=m:l_max;
lm=lm+1;
lmp=lmp+1;
aleg1(lm,n_theta) =plma(n_theta,lmp);
aleg2(lmp,n_theta)=2*pi*W(n_theta)*plma(n_theta,lmp);
end;
lmp=lmp+1;
aleg2(lmp,n_theta)   =2*pi*W(n_theta)*plma(n_theta,lmp);
end;
end;

% Now define what's useful for spectral to spatial transform

n_m_max=fix(m_max/minc)+1;
        lstart(1)=1;
        lstop(1)=l_max+2;
        lstrt(1)=1;
        lstp(1)=l_max+1;
        if ( mod(l_max,2)==0 )
           lmodd(1)=true;
        else
           lmodd(1)=false;
        end

        rim(1)=0.0;
        for mca=2:n_m_max
           mc=mca-1;
           lstart(mca)=lstop(mc)+1;
           lstop(mca)=lstart(mca)+l_max-mc*minc+1;
           lstrt(mca)=lstp(mc)+1;
           lstp(mca)=lstrt(mca)+l_max-mc*minc;
           rim(mca)=double(mc*minc);
           if ( mod(lstp(mca)-lstrt(mca),2)==0 )
              lmodd(mca)=true;
           else
              lmodd(mca)=false;
           end;
        end;
