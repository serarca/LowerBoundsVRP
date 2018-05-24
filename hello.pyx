def sum(x, y):
    print("Going to calculate sum of %d and %d..." % (x, y))
    return x + y

def a():
    x = 5
    b(x)

def b(x):
    y = 3
    total = sum(x, y)
    print("total = %d" % (total,))
