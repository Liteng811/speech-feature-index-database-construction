function f=py_method(filepath)
% str 输入wav音频文件的路径
% 
% 输出音频文件的语音识别结果
if count(py.sys.path,'') == 0 % 添加系统扫描路径
    insert(py.sys.path,int32(0),'');
end;
f = py.ASR.ASR(filepath); %
f = char(f);
return 
