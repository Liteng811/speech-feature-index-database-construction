    function [value10,loc10]=findpeaksreverse(value1,loc1)
    if numel(loc1)<=1
        loc10=value1;
        value10=loc1;
    elseif numel(loc1)>1
        ool=find(diff(loc1)<0, 1);  %����ӦС��0������ǿգ�����Ҫ�任
        oof=find(abs(floor(loc1)-loc1)>0,1);  %����ӦС��0������ǿգ�������С��������Ҫ�任
        oov=find(diff(value1)<0, 1);%��Ӧ��С��0���������û�У���ܿ��������,����У���˳���ǶԵģ����û���
        if ~isempty(ool) || ~isempty(oof)  %&& isempty(oov) 
            loc10=value1;
            value10=loc1;        
        elseif  ~isequal(round(loc1),loc1)     %�������ģ�һ����value ����loc
                loc10=value1;
                value10=loc1;
        elseif ~isempty(oov)
                value10=value1;
                loc10=loc1;
        else
                value10=value1;
                loc10=loc1;
        end
    end