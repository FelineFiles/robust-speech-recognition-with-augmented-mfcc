%% Main Program

%% Training
disp('####################  TRAINING   ###################');
tic
train(true)
toc


%% Testing
disp('####################  TESTING   ###################');
noise_levels = {'CLEAN', 'SNR10', 'SNR5'}
genders = {'male', 'female'}
tic
for k=1:numel(noise_levels)
    for l=1:numel(genders)
        test(noise_levels{k}, genders{l}, true)
    end
end
toc
