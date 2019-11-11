from scipy import linalg
import numpy as np
import sys
import inspect

A = None
B = None
C = None
X = None
Y = None
a = None
b = None
transA = False
transB = False

def main():
    var  = "ABCXYab"
    with open(sys.argv[1], "r") as f:
        lines = f.readlines()
    tmp = []
    m, n = 0, 0
    common = ''
    context = inspect.currentframe().f_back.f_locals
    for line in lines:
        if not line or line == "": continue
        if line[0] in var:
            common = line[0]
            # Set the tempory var to empty list and reset the corresponding
            # variable
            tmp = context[line[0]] = []
            words=line.split()
            if line[0] == 'a' or line[0] == 'b':
                context[line[0]] = float(words[1])
                continue
            # Add the transposition information
            if len(words[0]) == 2 and words[0][1] == "t":
                context["trans" + line[0]] = True
            m, n = map(int, words[1:])
            continue
        words=line.split()
        tmp.append(list(map(float, words)))
    for v in var:
        if type(context[v]) == list:
            context[v] = np.array(context[v])
            for i, elt in enumerate(context[v]):
                if type(elt) == list:
                    context[v][i] = np.array(elt)
    fname = sys.argv[2]
    func = getattr(linalg.blas, fname)
    if (fname == "dgemv"):
        print("transA ", transA)
        if transA:
            res = func(a, np.transpose(A), X, b, Y)
        else:
            res = func(a, A, X, b, Y)
        print("lengths : X : {}, Y : {}, res : {}".format(len(X), len(Y), len(res)))
    atol = 1e-9
    correct = np.allclose(res, C, rtol=0.0, atol=atol)
    print ("result of \"{}\" is correct : {}".format(fname, correct))
    for i in range(len(res)):
        if abs(res[i] - C[i]) > atol:
            print("First different value at {} : ({}, {})"\
                  .format(i, res[i], C[i]))
            return

main()
