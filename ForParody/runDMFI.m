% make a list of the graphic output files

!/bin/ls -1 Gfiles/  | grep Gt > ftr

% open the list

fnid=fopen('ftr')
incr=0

% main loop over the graphic files

while 0==0
incr=incr+1

[sfname,count]=fscanf(fnid,'%s1')
fname=['Gfiles/',sfname]

parodyload
if incr==1
init_prepare_mesh
prepare_leg
end

scene
prepare_mesh_reduced

makefig_1hem_v6

t1=text(-2,0,-1.5,sfname)
set(t1,'color','w')
set(t1,'fontsize',16)

fname=['Frames/frame',num2str(incr,'%0.3d')]
print (gcf,'-djpeg70',['Frames/frame',num2str(incr,'%0.3d')])

close(f)

end
