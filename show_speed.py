import multiprocessing as mp
import time
def compute(i):
    print("starting with number %s" % i)
    begin = time.time()
    if i % 2 == 0:
        print("finished")
        return(i*i,time.time()-begin)
    else:
        j = 0
        for s in range(0,100000000):
            j = j + i*i

        print("finished")
        return(j,time.time()-begin)


if __name__ == "__main__":
        p = mp.Pool(processes=3)
        results = p.map(compute,[1,2,1,2,1,2,1,2])
        for i in results:
            print(i)

