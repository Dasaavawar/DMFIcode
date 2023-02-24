% This routine performs inverse spherical harmonic transform over say flm2

s=flm2;

icn=0;
for icc=1:2:n_theta_max
ic1=icc+1;
icn=icn+1;
for n_m=1:n_m_max
lms=lstp(n_m);
s12=0.d0;
z12=0.d0;
for lm=lstrt(n_m):2:lms-1  
s12=s12+s(lm)*aleg1(lm,icn);
z12=z12+s(lm+1)*aleg1(lm+1,icn);
end;
if ( lmodd(n_m) )
     s12=s12+s(lms)*aleg1(lms,icn);
     end;
     sc(n_m,icc)=s12+z12;
     sc(n_m,ic1)=s12-z12;
     end;
     end;

% unscramble

sc2(1:n_m_max,1:n_theta_max/2)=sc(:,1:2:n_theta_max);
sc2(1:n_m_max,n_theta_max/2+1:n_theta_max)=sc(:,n_theta_max:-2:2);

% symetrise and normalize for Fourier

sc2(2:end,:)=sc2(2:end,:)/2;
sc2(n_m_max+1:n_phi_max/2+1,:)=0;
sc2(n_phi_max/2+2:n_phi_max,:)=conj(sc2(n_phi_max/2:-1:2,:));
sc2=ifft(sc2)*np;
