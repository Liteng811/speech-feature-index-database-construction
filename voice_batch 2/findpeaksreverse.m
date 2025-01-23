    function [value10,loc10]=findpeaksreverse(value1,loc1)
    if numel(loc1)<=1
        loc10=value1;
        value10=loc1;
    elseif numel(loc1)>1
        ool=find(diff(loc1)<0, 1);  %本不应小于0，如果非空，则需要变换
        oof=find(abs(floor(loc1)-loc1)>0,1);  %本不应小于0，如果非空，表明是小数，则需要变换
        oov=find(diff(value1)<0, 1);%本应有小于0的数，如果没有，则很可能是序号,如果有，则顺序是对的，不用换。
        if ~isempty(ool) || ~isempty(oof)  %&& isempty(oov) 
            loc10=value1;
            value10=loc1;        
        elseif  ~isequal(round(loc1),loc1)     %非整数的，一定是value 而非loc
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