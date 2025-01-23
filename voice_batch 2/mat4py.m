% pyversion C:\Users\wx\anaconda3\envs\wenet\python.exe
% %初此加载的时候可以写但是二次加载会报错
if count(py.sys.path,'') == 0 %空字符串里填写py文件的相对路径
    insert(py.sys.path,int32(0),'');
end;
t = py.ASR.ASR('mrher_short_.wav'); %
t = char(t);
disp(t)