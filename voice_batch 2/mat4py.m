% pyversion C:\Users\wx\anaconda3\envs\wenet\python.exe
% %���˼��ص�ʱ�����д���Ƕ��μ��ػᱨ��
if count(py.sys.path,'') == 0 %���ַ�������дpy�ļ������·��
    insert(py.sys.path,int32(0),'');
end;
t = py.ASR.ASR('mrher_short_.wav'); %
t = char(t);
disp(t)