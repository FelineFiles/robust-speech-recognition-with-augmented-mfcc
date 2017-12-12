function new_data=open_wavfile(filename)

tempfile='temp.wav';
command =['!ByteSwap ' filename ' ' tempfile];
eval(command);

fid=fopen(tempfile,'rb');
new_data=fread(fid,inf,'int16');
fclose(fid);

new_data=new_data/(max(new_data));

delete(tempfile);
