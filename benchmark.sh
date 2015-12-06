echo 'Benchmarking begins'
echo 'Case 1:'
time ./transport i benchmarks/1_base.txt | grep real
echo 'Case 2:'
time ./transport i benchmarks/2_hom.txt | grep real
echo 'Case 3:'
time ./transport i benchmarks/3_S2.txt | grep real
echo 'Case 4:'
time ./transport i benchmarks/4_tolx2.txt | grep real
echo 'Case 5:'
time ./transport i benchmarks/5_tolx5.txt | grep real
echo 'Case 6:'
time ./transport i benchmarks/6_tolx10.txt | grep real
echo 'Case 7:'
time ./transport i benchmarks/7_tol_2.txt | grep real
echo 'Case 8:'
time ./transport i benchmarks/8_tol_5.txt | grep real
echo 'Case 9:'
time ./transport i benchmarks/9_tol_10.txt | grep real
echo 'Benchmarking complete, now go manually figure out the times cause youre too lazy parse the output properly'
