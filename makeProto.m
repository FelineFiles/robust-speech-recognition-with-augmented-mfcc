% function makeProto(protoDir, hmmtype, feature, dimen, nstate, MaxAsync)
%---- 'protoDir' is the directory to which the prototype HMM should be outputted
%---- 'dimen' is the dimension(s) of the observation vectors
        %if the model is a single-stream model, dimen will be a single integer
        %if the model is a multi-stream model, dimen will be a vector of integers
%----- 'hmmtype' is the HMM type(IHMM=Independent, PHMM=product, CHMM=coupled)
function makeProto(protoDir, hmmtype, feature, dimen, nstate, MaxAsync)
doDebug = 0;
if doDebug
    hmmtype = 'PHMM';
    protoDir = 'dummy';
    dimen = [10 45];
    nstate = [3 3];
    MaxAsync = 2;
end
name = 'proto';
if isempty(whos('nstate'))
    nstate = 3 * ones(size(dimen));
end
if isempty(whos('MaxAsync'))
    MaxAsync = 1;
end

if length(nstate) == 1
    NoState = nstate;
elseif length(nstate) == 2
    statecomb = [];
    for i = 1: nstate(1)
        for j = 1: nstate(2)
            statecomb = [ statecomb; i j];
        end
    end
    idx = find(abs(statecomb(:,1)-statecomb(:,2))>MaxAsync);
    statecomb(idx,:) = [];
    NoState = length(statecomb);
else
    fprintf('[makeProto]ERROR: cannot handle multistream HMM with more than 2 streams states\n');
    return;
end
% observation probability parameters
mu     = zeros(NoState, sum(dimen));
sigma  =  ones(NoState,  sum(dimen));
gconst = zeros(NoState, length(dimen));;
% transition probabilities
switch hmmtype
case {'IHMM','MHMM'}
    transp = zeros(NoState+2);
    for i = 2: (NoState+1)
        transp(i, i)   = 0.6;
        transp(i, i+1) = 0.4;
    end
    transp(1, 2) = 1.0;
case {'PHMM' , 'CHMM'}
    transp = zeros(NoState+2);
    transp(1,2) = 1.0;
    for i = 2: (NoState+1)
        csc = statecomb(i-1, :);  % current state combination
        nsc = [csc; csc+[1 0]; csc+[0 1]; csc+[1 1]]; % next state combination
        nscidx = zeros(size(nsc,1),1);
        for j = 1: size(nsc,1) 
            t = find(statecomb(:,1) == nsc(j, 1) & statecomb(:,2) == nsc(j, 2) );
            if ~isempty(t)
                nscidx(j) = t;
            end
        end
        idx = find(nscidx ~= 0);
        nscidx = nscidx(idx) + 1;
        nsc    = nsc(idx, :);
        for j = 1: length(nscidx)
            if nscidx(j) == i
                transp(i, nscidx(j)) = 1 / length(nscidx) + 0.1;
            else
                transp(i, nscidx(j)) = (0.9 - 1/length(nscidx)) / (length(nscidx)-1);
            end
        end
    end
    transp(NoState+1, NoState+1) = 0.6;
    transp(NoState+1, NoState+2) = 0.4;
end
% output file
writeHMMdef([protoDir name], {name}, feature, dimen, {mu}, {sigma}, {transp}, {gconst}); 