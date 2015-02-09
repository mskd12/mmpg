#!/bin/bash

time=Results/mmpg_results_nash.txt
ExampleFolder=../MyExamples/MMPG

rm $time
echo -e "Input \t\t\t\t MMPG with Nash \t MMPG without Bribe \t MMPG with Bribe" > $time


for (( n=10; n <= 20; n++ ))
do
	for (( pl=5; pl <= 5; pl++ ))
	do
		for (( d=1; d <= 5; d++ ))
		do
			a=1;
			b=$((a+1))
			c=1
			#r=6
			example="randomgame-multi-$n-$pl-$d.pg"
			./pgsolver/bin/randomgame $n $c $a $b | ./src/mmpggenerator $pl > $ExampleFolder/$example
			#./pgsolver/bin/stratimprgen -pg switchbestexp $d | ./src/mmpggenerator > $ExampleFolder/$example
			#./pgsolver/bin/laddergame $d | ./src/mmpggenerator > $ExampleFolder/$example
			#./pgsolver/bin/steadygame $n 2 4 3 5 | ./src/mmpggenerator > $ExampleFolder/$example
			#n=$((n+1))
			mmpgNash="$(./src/mpg $ExampleFolder/${example} -b 2 -n $pl -o | grep FinalAnswer)"
			mmpgNormal="$(./src/mpg $ExampleFolder/${example} -b 0 -n $pl -o | grep FinalAnswer)"
			mmpgBribe="$(./src/mpg $ExampleFolder/${example} -b 1 -n $pl -o | grep FinalAnswer)"

			output="$example \t $mmpgNash \t $mmpgNormal \t $mmpgBribe"
			echo -e $output >> $time
		done
	done
done
