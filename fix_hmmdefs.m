function fix_hmmdefs(num_features,Root)

ModelRoot=[Root 'Models\'];
hmmdefs=[ModelRoot 'hmmdefs'];
sil_model=[ModelRoot 'sil'];

fid=fopen(hmmdefs);
ctr=1;
while(1)
    temp=fgets(fid);
    if(length(temp)==1 & temp==-1)
        break;
    else
        hmm{ctr}=temp;
        ctr=ctr+1;
    end;
end;
fclose(fid);

fid=fopen(sil_model);
ctr=1;
while(1)
    temp=fgets(fid);
    if(length(temp)==1 & temp==-1)
        break;
    else
        sil{ctr}=temp;
        ctr=ctr+1;
    end;
end;
fclose(fid);

for k=1:length(sil)
    if(length(sil{k})==10)
        if(sil{k}(1:9)=='<STATE> 8')
            mean_line=sil{k+2};
            var_line=sil{k+4};
            gconst_line=sil{k+5};
            break;
        end;
    end;
end;

fid=fopen(hmmdefs,'w');
for k=1:length(hmm)
    fprintf(fid,'%s',hmm{k});
end;
fprintf(fid,'%s\n',['~h "sp"']);
fprintf(fid,'%s\n',['<BEGINHMM>']);
fprintf(fid,'%s\n',['<NUMSTATES> 5']);
for state=2:4
    fprintf(fid,'%s\n',['<STATE> ' num2str(state)]); 
    fprintf(fid,'%s\n',['<MEAN> ' num2str(num_features)]);
    fprintf(fid,'%s',mean_line);
    fprintf(fid,'%s\n',['<VARIANCE> ' num2str(num_features)]);
    fprintf(fid,'%s',var_line);
end;
fprintf(fid,'%s\n',['<TRANSP> ' num2str(5)]);
fprintf(fid,'%s\n',[num2str(0) ' ' num2str(1) ' ' num2str(0) ' ' num2str(0) ' ' num2str(0)]);
fprintf(fid,'%s\n',[num2str(0) ' ' num2str(.5) ' ' num2str(.5) ' ' num2str(0) ' ' num2str(0)]);
fprintf(fid,'%s\n',[num2str(0) ' ' num2str(0) ' ' num2str(.5) ' ' num2str(.5) ' ' num2str(0)]);
fprintf(fid,'%s\n',[num2str(0) ' ' num2str(0) ' ' num2str(0) ' ' num2str(.5) ' ' num2str(.5)]);
fprintf(fid,'%s\n',[num2str(0) ' ' num2str(0) ' ' num2str(0) ' ' num2str(0) ' ' num2str(0)]);
fprintf(fid,'%s\n',['<ENDHMM>']);
fclose(fid);
    
    

