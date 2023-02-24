% This routine performs space to spectral over say Vr

anlc=fft(sc2)/4/pi/np;

%% scrambling

for n_theta=1:n_theta_max/2
anlc1(:,2*n_theta-1)=anlc(:,n_theta);
anlc1(:,2*n_theta)=anlc(:,n_theta_max-n_theta+1);
end;

n_theta_hs=0;
for n_theta_n=1:2:n_theta_max;
n_theta_s=n_theta_n+1;   
n_theta_hs=n_theta_hs+1;

for mca=1:n_m_max   
a1plus( mca,n_theta_hs)=anlc1(mca,n_theta_n) + anlc1(mca,n_theta_s);
a1minus(mca,n_theta_hs)=anlc1(mca,n_theta_n) - anlc1(mca,n_theta_s);
end;
end;

for mca=1:n_m_max

lms=lstop(mca);
lm1=lms-1;
ap_1=a1plus(mca,1);
ap_2=a1plus(mca,2);
am_1=a1minus(mca,1);
am_2=a1minus(mca,2);
for lmp=lstart(mca):2:lms-1
lm1=lmp+1;
flm1(lmp)=ap_1*aleg2(lmp,1)+ap_2*aleg2(lmp,2);
flm1(lm1)=am_1*aleg2(lm1,1)+am_2*aleg2(lm1,2);
end;

if ( lm1<lms )
     flm1(lms)=ap_1*aleg2(lms,1)+ap_2*aleg2(lms,2);
     end;

     end;   

n_theta_1=1;

for n_theta_rel_1=3:2:n_theta_max/2
n_theta_rel_2=n_theta_rel_1+1;
n_theta_1=n_theta_1+2;
n_theta_2=n_theta_1+1;

for mca=1:n_m_max;
lms=lstop(mca);
lm1=lms-1;
ap_1=a1plus(mca,n_theta_rel_1);
ap_2=a1plus(mca,n_theta_rel_2);
am_1=a1minus(mca,n_theta_rel_1);
am_2=a1minus(mca,n_theta_rel_2);
for lmp=lstart(mca):2:lms-1;
lm1=lmp+1;
flm1(lmp)=flm1(lmp) + ap_1*aleg2(lmp,n_theta_1) + ap_2*aleg2(lmp,n_theta_2);
flm1(lm1)=flm1(lm1) + am_1*aleg2(lm1,n_theta_1) + am_2*aleg2(lm1,n_theta_2);
end;

              if ( lm1<lms )
flm1(lms)=flm1(lms) + ap_1*aleg2(lms,n_theta_1) + ap_2*aleg2(lms,n_theta_2);
end;

end;
end;

clear flm2
for mca=1:n_m_max;
flm2(lstrt(mca):lstp(mca))=flm1(lstart(mca):lstop(mca)-1);
end
