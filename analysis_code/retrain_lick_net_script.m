%retrain_lick_net_script
%animals = {'bb8','bb9','bb10','bb11','bb12','bb13','bb14','bb15','bb16','bb17','bb18','bb19'};
%animals = {'bb10','bb11','bb12','bb13','bb14','bb15','bb16','bb17','bb18','bb19','mt1','mt2'};
%dates1 = {'190629','190630','190701','190702','190703','190704','190705','190808','190809','190810','190811','190812','190813','190814'};
%dates2 = {'190816','190817','190818','190819','190820','190821','190822','190823','190824','190825','190826','190827','190828'};
%dates3 = {'191027','191028','191029','191030','191031','191101','191102','191103','191104','191105','191106','191107','191108','191109'};
%dates4 = {'200124','200125','200126','200127','200128','200129','200130','200131','200201','200202','200203','200204','200205','200206','200207'};
%dates5 = {'200209','200210','200211','200212','200213','200214','200215','200216','200217','200218','200219','200220','200221','200222'};
%allDates = [dates1 dates2 dates3 dates4 dates5];
dateInfo = load('data/dateInfo.mat'); dateInfo=dateInfo.dateInfo;

animals = {dateInfo.animal};
allDates = {};
count = 0;
for i=1:length(dateInfo)
    for j=1:length(dateInfo(i).dates)
        hasDate = false;
        for k=1:length(allDates)
            if (strcmp(dateInfo(i).dates{j},allDates{k}))
                hasDate = true;
            end
        end
        if (~hasDate)
            count = count + 1;
            allDates{count} = dateInfo(i).dates{j};
        end
    end
end
[lickNet] = retrainLickNet_subsample('2bottle',allDates,animals,10000);