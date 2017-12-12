% writeHMMdef.m, by J. Xue, -20050915-
%  generate HMM definition files
%function writeHMMdef(filename, names, dimen, mu, sigma, transp, gconst)
%  Only single mixture case is implemented, gconst is optional, other inputs are required

function writeHMMdef(filename, names, feature, dimen, mu, sigma, transp, gconst )
doDebug = 0;
if doDebug
    filename = 'dummy';
    names = {'proto'};
    dimen = [3 3];
    state = 5;
    mu    = {zeros(state, sum(dimen))};
    sigma = {zeros(state, sum(dimen))};
    gconst= {zeros(state, length(dimen))};
    transp= {zeros(state+2, state+2)};
end

NoState = size(transp{1},1);
if length(dimen)>1
    sweight = ones(length(dimen),1);
end
col = cell(length(dimen),1);
for i = 1: length(col)
    if i == 1
        col{i} = [1:dimen(1)];
    else
        col{i} = [sum(dimen(1:(i-1)))+1 : sum(dimen(1:i))];
    end
end
fid = fopen(filename, 'w');
fprintf(fid, '~o\n');
fprintf(fid, '<VECSIZE> %d <NULLD><%s>\n', sum(dimen), feature);
if length(dimen) > 1
    fprintf(fid, '<STREAMINFO> %d %s\n', length(dimen), num2str(dimen(:)'));
end
for i = 1: length(names)
    fprintf(fid, '~h "%s"\n', names{i});
    fprintf(fid, '<BEGINHMM>\n<NUMSTATES> %d\n', NoState);
    for j = 2: (NoState-1)
        fprintf(fid, '<STATE> %d\n', j);
        if length(dimen) > 1
            fprintf(fid, '<SWEIGHTS> %d %s\n', length(sweight), num2str(sweight(:)'));
        end
        for k = 1: length(dimen)
            if length(dimen) > 1
                fprintf(fid, '<STREAM> %d\n', k);
            end
            fprintf(fid, '<MEAN> %d\n', dimen(k));
            for l = 1: length(col{k})
                fprintf(fid, '%.6f ', mu{i}(j-1, col{k}(l)) );
            end
            fprintf(fid, '\n');
            fprintf(fid, '<VARIANCE> %d\n', dimen(k));
            for l = 1: length(col{k})
                fprintf(fid, '%.6f ', sigma{i}(j-1, col{k}(l)) );
            end
            fprintf(fid, '\n');
            if ~isempty(whos('gconst'))
                fprintf(fid, '<GCONST> %d\n', gconst{i}(j-1,k) );
            end
        end
    end
    fprintf(fid, '<TRANSP> %d\n', NoState);
    for j = 1: NoState
        for k = 1: NoState
            fprintf(fid, '%.6f ', transp{i}(j, k));
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, '<ENDHMM>\n');
end
fclose(fid);