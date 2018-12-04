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

    void ASKI(const string &prompt, int &iinput);
    void ASKR(const string &prompt, double &rinput);
    void ASKL(const string &prompt, bool &linput);
    void ASKS(const string &prompt, string &input);
    void ASKC(const string &prompt, string &command, string &cargs);

    void LC2UC(string &input);

    void READI(int n, ivec &ivar, bool &error);
    void READR(int n, vec &var, bool &error);

    void GETINT(string &input, ivec &a, int &n, bool &error);
    void GETFLT(string &input, vec &a, int &n, bool &error);
}

#endif //XROTOR_NOGRAPHICS_CPP_USERIO_H
