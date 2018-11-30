//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_USERIO_H
#define XROTOR_NOGRAPHICS_CPP_USERIO_H

#include <vector>   // std::vector
#include <string>   // std::getline, std::string, std::stoi, std::stod
using std::string;
using std::vector;

typedef vector<int> ivec;
typedef vector<double> vec;

namespace userio {

    void aski(const string &prompt, int &iinput);
    void askr(const string &prompt, double &rinput);
    void askl(const string &prompt, bool &linput);
    void asks(const string &prompt, string &input);
    void askc(const string &prompt, string &command, string &cargs);

    void lc2uc(string &input);

    void readi(int n, ivec &ivar, bool &error);
    void readr(int n, vec &var, bool &error);

    void getint(string &input, ivec &a, int &n, bool &error);
    void getflt(string &input, vec &a, int &n, bool &error);
}

#endif //XROTOR_NOGRAPHICS_CPP_USERIO_H
