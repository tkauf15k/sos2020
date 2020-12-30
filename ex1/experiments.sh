outbasedir="./output/simple"
mkdir -p $outbasedir
algorithm=ga
timelimit=300
localsearch=False
instances=(
10
20
30
40
50
60
70
80
90
100 
)
# 125
# 150
# 175
# 200
# )

for i in "${instances[@]}";
do
	for j in {0..9}; 
	do
		outfile="${outbasedir}/resultdd_${i}_${j}.txt"
		echo $outfile
		(python main_evaluator.py --instance $i --time $timelimit --output $outfile $algorithm -l $localsearch)&
		pids[${i}]=$!
	done
	echo "waiting..."
	for pid in ${pids[*]}; do
		wait $pid
	done
	echo "done waiting...";
done