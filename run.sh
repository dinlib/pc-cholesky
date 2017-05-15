#! /bin/bash

echo "Polybench Cholesky. Sequential x OpenMP parallel x Pthread parallel"

# CC=/usr/local/bin/gcc-6
CC=gcc
CFLAGS="-std=gnu99 -lm"
DATASET="-DLARGE_DATASET"
nthreads=32

# SEQUENTIAL CHOLESKY
# echo "SEQUENTIAL"
$CC $CFLAGS -I utilities utilities/polybench.c $1.c $DATASET -o $1.out
file="./data/$1_sequential.data"
if [ -f "$file" ]
then
	rm $file
fi
for i in {1..11}
do
	 echo "SEQUENTIAL - EXECUTION $i"
   if [[ $i == 1 ]]; then
     ./$1.out > /dev/null
   else
    # ./$1.out $MSIZE 2 >> testfile.data
     ./$1.out >> ./data/$1_sequential.data 2>&1
  fi
done
cat ./data/$1_sequential.data | sort > sort.data
sed '1d; $d' sort.data > final.data
avgs=$(awk '{s+=$1}END{print s/NR}' RS="\n" final.data)
rm final.data sort.data
# ./cholesky.out

# OMP CHOLESKY
$CC $CFLAGS -fopenmp -I utilities utilities/polybench.c $1_omp.c $DATASET -o $1_omp.out
avgomp=(0 0 0 0 0)
for t in 2 4 8 16 32
do
  # echo "OMP FOR $t THREADS"
  file="./data/$1_omp_$t.data"
  if [ -f "$file" ]
  then
  	rm $file
  fi
  for i in {1..11}
  do
		 echo "OMP FOR $t THREADS - EXECUTION $i"
     if [[ $i == 1 ]]; then
       ./$1_omp.out $t > /dev/null
     else
      # ./$1.out $MSIZE 2 >> testfile.data
       ./$1_omp.out $t >> ./data/$1_omp_$t.data 2>&1
    fi
  done
  cat ./data/$1_omp_$t.data | sort > sort.data
  sed '1d; $d' sort.data > final.data
  avgomp[$((t-1))]=$(awk '{s+=$1}END{print s/NR}' RS="\n" final.data)
  echo ${avgomp[$((t-1))]}
  rm final.data sort.data
done

# PTHREAD CHOLESKY
$CC $CFLAGS -pthread -I utilities utilities/polybench.c $1_pthread.c $DATASET -o $1_pthread.out
avgp=(0 0 0 0 0)
for t in 2 4 8 16 32
do
	file="./data/$1_pthread_$t.data"
	if [ -f "$file" ]
	then
		rm $file
	fi
	# echo "PTHREAD FOR $t THREADS"
	for i in {1..11}
	do
		echo "PTHREAD FOR $t THREADS - EXECUTION $i"
		if [[ $i == 1 ]]; then
			./$1_pthread.out $t > /dev/null
		else
			# ./$1.out $MSIZE 2 >> testfile.data
			./$1_pthread.out $t >> ./data/$1_pthread_$t.data 2>&1
		fi
	done
	cat ./data/$1_pthread_$t.data | sort > sort.data
	sed '1d; $d' sort.data > final.data
	avgp[$((t-1))]=$(awk '{s+=$1}END{print s/NR}' RS="\n" final.data)
	echo ${avgp[$((t-1))]}
	rm final.data sort.data
done

#SAVING FILE FOR PLOTING IN GNUPLOT
file="./data/$1_speedup.data"
if [ -f $file ] ; then
    rm $file
fi
for ((t=0; t<nthreads; t+2)); do
	if [[ $t == 0 ]]; then
		echo "1" "1"  "1" >> ./data/$1_speedup.data
	else
		pthread=$(echo "$avgs / ${avgp[t]}" | bc -l)
		omp=$(echo "$avgs / ${avgomp[t]}" | bc -l)
		echo $((t+1)) $omp $pthread >> ./data/$1_speedup.data
 fi
done
