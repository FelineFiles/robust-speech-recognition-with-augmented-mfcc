function test(noise_level, gender, do_disp)
    if(do_disp)  
        disp('=================================================');
        disp(['==================  ', noise_level, gender, ' =================']);
        disp('=================================================');
    end;
    Fs=8000;

    Root='.\';

    TrainDataRoot=[Root 'database\train\'];
    TestDataRoot=[Root 'database\test\', noise_level, '_', gender, '\'];

    TrainFeatureRoot=[Root 'features\train\'];
    TestFeatureRoot=[Root 'features\test\'];

    ModelRoot=[Root 'models\'];
    ScriptRoot=[Root 'scripts\'];
    wdnet_file      = [ScriptRoot 'wdnet_file.txt'];

    ConfigFile=[ScriptRoot 'config_file.cfg'];
    ModelList=[ScriptRoot 'Model_List.txt'];
    DictFile=[ScriptRoot 'Dictionary.txt'];
    WordList=[ScriptRoot 'Word_List.txt'];
    WordList2=[ScriptRoot 'Word_List2.txt'];
    WordListSP=[Root 'Scripts\WordListSP.txt'];
    MLF_Results=[ScriptRoot 'MLF_Results.mlf'];
    TrainWordMLF=[Root 'Scripts\TrainWordMLF.mlf'];
    TrainWordMLF2=[Root 'Scripts\TrainWordMLF2.mlf'];
    TestWordMLF=[Root 'Scripts\TestWordMLF.mlf'];
    TrainFeatureScript=[Root 'Scripts\TrainFeatureScript.txt'];
    TestFeatureScript=[Root 'Scripts\TestFeatureScript.txt'];
    TestScript=[Root 'Scripts\TestScript.txt'];
    MixScript1=[Root 'Scripts\HED1.txt'];
    MixScript2=[Root 'Scripts\HED2.txt'];
    WdnetFile=[Root 'Scripts\WDNet.txt'];
    MLFResults=[Root 'Scripts\MLFResults.mlf'];
    hmmdefs=[ModelRoot 'hmmdefs'];

    NUM_STATES=16;
    NUM_HEREST_1=3;
    NUM_HEREST_2=6;

    % =======================================================================================
    % =============== Testing Feature Extraction ============================================
    % =======================================================================================
    testfiles=dir(TestDataRoot);
    testfiles=testfiles(3:end);
    features=dir(TestFeatureRoot);
    for n=3:length(features)
        delete([TestFeatureRoot '\' features(n).name]);  
    end;
    if(do_disp)

        disp('>>> Performing testing feature extraction');

    end;     
    numfiles=length(testfiles);
    for num=1:numfiles
        if(do_disp & mod(num,200)==0)
            disp([num2str(ceil(100*num/numfiles)) '% done...']);
        end;
        file_name=char(testfiles(num).name);
        wavFile=[TestDataRoot file_name];
        load(wavFile);
        data = new_data;
        feature=extract_features(data,Fs, 'train');
        feature_file=[TestFeatureRoot file_name(1:end-2) 'mfc'];
        writehtk(feature_file,feature,1/120,9);
    end;
    fid=fopen(TestFeatureScript,'w');
    for i=1:numfiles
        fprintf(fid, '%s\n',[TestFeatureRoot testfiles(i).name(1:end-2) 'mfc']);
    end
    fclose(fid);

    if(do_disp)  

        disp('>>> Testing HMMs');

    end;
        
    disp(['Creating MLF file...']);  
    feature_files=char(textread(TestFeatureScript,'%s'));
    fid1=fopen(TestWordMLF,'w');
    fprintf(fid1,'%s\n','#!MLF!#');
    for k=1:size(feature_files,1)
        dashes=find(feature_files(k,:)=='-');
        dots=find(feature_files(k,:)=='.');
        slashes=find(feature_files(k,:)=='\');
        underscores=find(feature_files(k,:)=='_');
        for s=1:length(slashes)
            feature_files(k,slashes(s))='/';
        end;  
        fprintf(fid1,'%s\n',['"' feature_files(k,1:dots(end)) 'lab"']);
        fprintf(fid1,'%s\n','sil');
        words=feature_files(k,underscores(end)+1:dots(end)-2);
        for w=1:length(words)
            number=find_number(words(w));
            fprintf(fid1,'%s\n',number);
        end;
        fprintf(fid1,'%s\n','sil');
        fprintf(fid1,'%s\n','.');       
    end;
    fclose(fid1);
    cmd=['!HBuild'...
        ' -s sil sil'...
        ' ' WordList...
        ' ' WdnetFile...
        ];
    eval(cmd);
    disp('HVite.');
    cmd = ['!HVite' ...
        ' -C ' ConfigFile...
        ' -H ' hmmdefs...
        ' -i ' MLFResults...
        ' -I ' TestWordMLF...
        ' -w ' WdnetFile...
        ' -p -20.0 '...
        ' -S ' TestFeatureScript ...
        ' ' DictFile ...
        ' ' WordList ...
        ];
    eval(cmd);        
    disp('HREsults.');
    cmd = ['! HResults '...
        '   -e  "???" sil'...
        ' -I ' TestWordMLF...
        ' -p '...
        ' ' WordList...
        ' ' MLFResults...
        ];
    eval(cmd);

end