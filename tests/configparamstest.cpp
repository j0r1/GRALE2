#include "configurationparameters.h"
#include <serut/vectorserializer.h>
#include <iostream>
#include <fstream>

using namespace grale;
using namespace std;
using namespace serut;

void writeBuf(const string fname, const vector<uint8_t> &data)
{
    ofstream f(fname, ios::out | ios::binary);
    f.write((char *)data.data(), data.size());
}

int main(int argc, char const *argv[])
{
    int errCount = 0;
    auto test = [&errCount](const string &name, auto param)
    {
        cerr << "Test " << name << endl;
        VectorSerializer ser;

        {
            cerr << "  Out:" << endl;
            TypedParameter t { param };
            t.dump("    ");

            if (!t.write(ser))
            {
                cerr << "ERROR: Couldn't write data: " << t.getErrorString() << endl;
                errCount++;
                return;
            }
        }

        VectorSerializer ser2;

        {
            TypedParameter t;
            if (!t.read(ser))
            {
                cerr << "ERROR: Couldn't read data: " << t.getErrorString() << endl;
                errCount++;
                return;
            }
            cerr << "  In:" << endl;
            t.dump("    ");

            if (!t.write(ser2))
            {
                cerr << "ERROR: Couldn't write second data: " << t.getErrorString() << endl;
                errCount++;
                return;
            }
        }

        if (ser2.getBuffer() != ser.getBuffer())
        {
            cerr << "ERROR: Buffers don't match" << endl;
            writeBuf(name + "buf1.dat", ser.getBuffer());
            writeBuf(name + "buf2.dat", ser2.getBuffer());
            errCount++;
            return;
        }
    };

    test("bool", true);
    test("int", 123);
    test("real", 45.6);
    test("string", string("Hello world!"));

    test("bool", vector<bool>{true, false});
    test("int", vector<int>{123, 456});
    test("real", vector<double>{45.6, 133.7});
    test("string", vector<string>{"Hey", "there"});


    cerr << "Detected " << errCount << " errors" << endl;    
    return 0;
}
