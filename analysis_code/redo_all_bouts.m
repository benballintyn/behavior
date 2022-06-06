%redo_all_bouts
animals = {'bb8','bb9','bb10','bb11','bb12','bb13'};
dates1 = {'190629','190630','190701','190702','190703','190704','190705','190808','190809','190810','190811','190812','190813','190814'};
dates2 = {'190816','190817','190818','190819','190820','190821','190822','190823','190824','190825','190826','190827','190828'};
dates{1} = dates1; dates{2}=dates1; dates{3} = dates2; dates{4} = dates2; dates{5} = dates2; dates{6} = dates2;

for i=1:length(animals)
    for j=1:length(dates{i})
        disp(['Date: ' dates{i}{j} ' animal: ' animals{i}])
        redoBouts('analyzed_data',dates{i}{j},animals{i});
    end
end