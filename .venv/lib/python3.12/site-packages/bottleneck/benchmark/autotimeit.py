import timeit


def autotimeit(stmt, setup="pass", repeat=3, mintime=0.2):
    timer = timeit.Timer(stmt, setup)
    number, time1 = autoscaler(timer, mintime)
    time2 = timer.repeat(repeat=repeat - 1, number=number)
    return min(time2 + [time1]) / number


def autoscaler(timer, mintime):
    number = 1
    for i in range(12):
        time = timer.timeit(number)
        if time > mintime:
            return number, time
        number *= 10
    raise RuntimeError("function is too fast to test")
