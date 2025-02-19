cmake -DCMAKE_C_COMPILER=/home/lixinyan/gcc-10/bin/gcc -DCMAKE_CXX_COMPILER=/home/lixinyan/gcc-10/bin/g++ ..
make -j36
cd bin

cd output_temp
rm -r ./*
cd ..

cd simulate_data
rm -r ./*
cd ..

cd simulate_data_noOne
rm -r ./*
cd ..

rm -rf error_LM_Smoothing.txt
rm -rf error_LM.txt
rm -rf error_Smoothing_LM.txt
rm -rf error_Smoothing.txt
rm -rf memory.txt
rm -rf  error_inner_py.txt
rm -rf error_begin.txt

gcc run_testmemory.c -o run_testmemory

for((i=0;i<1;i++));do

    echo "--------------------------------- "$[i]" -------------------------------------------"
    for((j=0;j<1;j++));do
    ./run_testmemory $[j] $[i]
    done
    for((j=4;j<5;j++));do
    ./run_testmemory $[j] $[i]
    done
    for((j=6;j<8;j++));do
    ./run_testmemory $[j] $[i]
    done
   ######python3 BBb_L1_interiorpoint.py $[i]
done

cd ..


