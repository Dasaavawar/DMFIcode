
lm=0;
for mca=0:n_m_max-1;
for l=azsym*mca:l_max;
lm=lm+1;

flm2(lm)=flm2(lm)*0.536^(l+2);
end;
end;

