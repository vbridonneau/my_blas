double my_ddot(const int N, const double *X, const int incX, const double *Y, const int incY) {
    double res = 0.0;
    int xcpt = 0, ycpt = 0;
    for (int elt = 0; elt < N; elt++, xcpt += incX, ycpt += incY) {
        res += X[xcpt] * Y[ycpt];
    }
    return res;
}