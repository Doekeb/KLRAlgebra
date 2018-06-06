import threading
from threading import Thread
import time, random, Queue

def doubler(number):
    print(threading.currentThread().getName())
    print(number * 2)
    print('\n')

jack = 0
jill = 100

def inc():
    global counter
    counter += 1

def dec():
    global counter
    counter -= 1

def give_to_jack():
    global jack, jill
    time.sleep(random.random())
    if jill > 0:
        jill -= 1
        jack += 1

def give_to_jill():
    global jack, jill
    time.sleep(random.random())
    if jack > 0:
        jack -= 1
        jill += 1

def funnel_to_jack():
    while jill > 0:
        give_to_jack()

def funnel_to_jill():
    while jack > 0:
        give_to_jill

# if __name__ == '__main__':
#     for _ in range(10):
#         # time.sleep(random.random())
#         threading.Thread(target=funnel_to_jack).start()
#         threading.Thread(target=funnel_to_jill).start()





def summ(L):
    if L == []:
        return 0
    elif len(L) == 1:
        return L[0]
    else:
        time.sleep(1)
        return summ(L[:len(L)/2]) + summ(L[len(L)/2:])

def sumt(L):
    if L == []:
        return 0
    elif len(L) == 1:
        return L[0]
    else:
        time.sleep(1)
        q = Queue.Queue()
        left = Thread(target=lambda q, M: q.put(sumt(M)), args=(q, L[:len(L)/2]))
        right = Thread(target=lambda q, M: q.put(sumt(M)), args=(q, L[len(L)/2:]))
        left.start()
        right.start()
        left.join()
        right.join()
        return q.get() + q.get()


if __name__ == '__main__':
    L = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]

    t0 = time.time()
    print summ(L), time.time()-t0

    t1 = time.time()
    print sumt(L), time.time()-t1



# from threading import Thread
# from multiprocessing import Pool
# import time
#
# def show(x):
#     print x
#
# def show_n(x, n):
#     for _ in range(n):
#         time.sleep(.1)
#         print x
#
# def show_50(x):
#     show_n(x, 50)
#
# def show_forever(x):
#     while True:
#         print x
#
# def re(x):
#     return x
#
# counter = 0
#
# def increment_n(counter,n):
#     for _ in range(n):
#         time.sleep(.1)
#         counter += 1
#         print counter
#
# def decrement_n(counter,n):
#     for _ in range(n):
#         time.sleep(.1)
#         counter -= 1
#         print counter
#
# pool = Pool(2)
# thing = pool.map(show_50, ["hello", "goodbye"])
# thing

# hello = Thread(target=increment_n, args=(counter,50))
# goodbye = Thread(target=decrement_n, args=(counter,50))
#
# hello.start()
# goodbye.start()
