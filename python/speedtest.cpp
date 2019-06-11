#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include <iostream>

using namespace std;

int main(){
    constexpr size_t bufferSize = 10000;
    string str;
    auto buffer = (new char[bufferSize]);
    stringbuf ss;
    ifstream bigFile("data.pos_1.xyz");
    while (bigFile.good() and !bigFile.eof()){
        //bigFile.getline(buffer,bufferSize);
        //getline(bigFile,str);
        bigFile.get(ss);
        ss >> str;
        cout << str;
        //bigFile.read(buffer.get(), bufferSize);
        //bigFile.readsome(buffer.get(), bufferSize);
        //string tmps(buffer.get(),buffer.get()+bufferSize);
        //cout << tmps << "*********";
    }
}
/*
constexpr size_t bufferSize = 10 * 10;
unique_ptr<char[]> buffer(new char[bufferSize]);
while (bigFile)
{
    bigFile.read(buffer.get(), bufferSize);
    // process data in buffer
    string tmps(buffer.get(),buffer.get()+bufferSize);
    cout << tmps << "*********";
    }
*/
