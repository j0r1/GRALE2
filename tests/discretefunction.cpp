#include "discretefunction.h"
#include <iostream>

using namespace grale;
using namespace std;
using namespace errut;

void checkFalse(bool_t r, const string &desc)
{
    if (r)
        cout << "ERROR: expected failure for " << desc << ": " << r.getErrorString() << endl;
    else
        cout << "OK: got expected failure for " << desc << endl;
}

void checkTrue(bool_t r, const string &desc)
{
    if (!r)
        cout << "ERROR: expected success for " << desc << ": " << r.getErrorString() << endl;
    else
        cout << "OK: got expected success for " << desc << endl;
}

template <class T>
void checkValue(DiscreteFunction<T> &f, T x, T expectedY)
{
    T y = f(x);
    if (y != expectedY)
        cout << "ERROR: expected y=" << expectedY << " for x=" << x << " but found y=" << y 
             << "(difference is " << (y-expectedY) << ")" << endl;
    else
        cout << "OK: found y=" << y << " for x=" << x << endl;
}

template<class T>
void test()
{
    DiscreteFunction<T> f;

    auto r = f.init(1, 1, { (T)0, (T)0 });
    checkFalse(r, "same x0 and x1");

    r = f.init(1, 0, { (T)0, (T)0 });
    checkFalse(r, "x0 larger than x1");

    r = f.init(1, 0, { });
    checkFalse(r, "empty vector");

    r = f.init(1, 2, { (T)1 });
    checkFalse(r, "vector of length one");

    r = f.init(0, 1, { (T)1, (T)2 });
    checkTrue(r, "vector of length 2");
    checkValue(f, (T)0.75, (T)1.75);
    checkValue(f, (T)(-0.1), (T)1);
    checkValue(f, (T)(2.1), (T)2);

    r = f.init(10, 14, { (T)1, (T)2, (T)0, (T)4, (T)3 });
    checkTrue(r, "vector of length 5");
    checkValue(f, (T)9, (T)1);
    checkValue(f, (T)10.5, (T)1.5);
    checkValue(f, (T)11.75, (T)0.5);
    checkValue(f, (T)12.0, (T)0);
    checkValue(f, (T)13.25, (T)3.75);
    checkValue(f, (T)14, (T)3);
}

int main(int argc, char const *argv[])
{
    test<float>();
    test<double>();
    return 0;
}
