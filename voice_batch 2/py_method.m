function f=py_method(filepath)
% str ����wav��Ƶ�ļ���·��
% 
% �����Ƶ�ļ�������ʶ����
if count(py.sys.path,'') == 0 % ���ϵͳɨ��·��
    insert(py.sys.path,int32(0),'');
end;
f = py.ASR.ASR(filepath); %
f = char(f);
return 
