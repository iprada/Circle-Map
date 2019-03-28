

def generator():
    my_Dict = {}
    for i in range(0,100):
        print("generating %s" % i)
        my_Dict[i] = i
        yield

print(generator())


for j in generator():
    for i in range(0,10000000):
        a = 0
    print(j)