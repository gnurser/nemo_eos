thread_speed(){
    neos=$1
    nthreads=$2
    short_subroutine=$3
    subroutine=eos_${short_subroutine:=insitu4}
    echo neos=${neos}, nthreads=${nthreads}, subroutine=${subroutine}
    python -m timeit --setup "from set_test import set_test_values; fillvalue,mask,T4,S4,depth4 = set_test_values($length); from nemo_eos.nemo_rho import eos; eos.eos_init($neos); eos.set_eos_threads($nthreads)" "eos.${subroutine}(T4,S4,depth4)"
}
thread_speed_m(){
    neos=$1
    nthreads=$2
    short_subroutine=$3
    subroutine=eos_${short_subroutine:=insitu4_m}
    echo neos=${neos}, nthreads=${nthreads}, subroutine=${subroutine}
    python -m timeit --setup "from set_test import set_test_values; fillvalue,mask,T4,S4,depth4 = set_test_values($length); from nemo_eos.nemo_rho import eos; eos.eos_init($neos); eos.set_eos_threads($nthreads)" "eos.${subroutine}(fillvalue, mask, T4,S4,depth4)"
}
thread_speed_insitu04(){
    neos=$1
    nthreads=$2
    short_subroutine=insitu04
    subroutine=eos_${short_subroutine}
    echo neos=${neos}, nthreads=${nthreads}, subroutine=${subroutine}
    python -m timeit --setup "from set_test import set_test_values; fillvalue,mask,T4,S4,depth4 = set_test_values($length); from nemo_eos.nemo_rho import eos; eos.eos_init($neos); eos.set_eos_threads($nthreads)" "eos.${subroutine}(2.5,35.0,1.,depth4)"
}
length=$1
length=${length:=100000}

echo old EOS
thread_speed -1 1 insitu4
thread_speed -1 4 insitu4
thread_speed -1 8 insitu4
thread_speed_m -1 1 insitu4_m
thread_speed_m -1 4 insitu4_m

thread_speed_insitu04 -1 1 insitu04
thread_speed_insitu04 -1 4 insitu04

thread_speed -1 1 rab4
thread_speed -1 4 rab4
thread_speed_m -1 1 rab4_m
thread_speed_m -1 4 rab4_m

echo new EOS
thread_speed 2 1 insitu4
thread_speed 2 4 insitu4
thread_speed 2 8 insitu4
thread_speed_m 2 1 insitu4_m
thread_speed_m 2 4 insitu4_m

thread_speed_insitu04 2 1 insitu04
thread_speed_insitu04 2 4 insitu04

thread_speed 2 1 rab4
thread_speed 2 4 rab4
thread_speed_m 2 1 rab4_m
thread_speed_m 2 4 rab4_m





#     import timeit

#     print(timeit.timeit("fillvalue,mask,T4,S4,depth4 = set_test_values(100000)", setup="from set_test import set_test_values"))
     # print(timeit.timeit("eos.eos_insitu4(T4,S4,depth4)", setup="from nemo_eos.nemo_rho import eos; eos.eos_init(-1); eos.set_eos_threads(1); from set_test import set_test_values; fillvalue,mask,T4,S4,depth4 = set_test_values(100000)"))
     # timeit rho = eos.eos_insitu4_m(fillvalue, mask,T4,S4,depth4)
     # timeit alpha, beta = eos.eos_rab4(T4,S4,depth4)
     # timeit alpha, beta = eos.eos_rab4_m(fillvalue, mask,T4,S4,depth4)
     # eos.set_eos_threads(4)
     # timeit rho = eos.eos_insitu4(T4,S4,depth4)
     # timeit rho = eos.eos_insitu4_m(fillvalue, mask,T4,S4,depth4)
     # timeit alpha, beta = eos.eos_rab4(T4,S4,depth4)
     # timeit alpha, beta = eos.eos_rab4_m(fillvalue, mask,T4,S4,depth4)
