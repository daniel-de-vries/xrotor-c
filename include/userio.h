//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_USERIO_H
#define XROTOR_NOGRAPHICS_CPP_USERIO_H

#include <string>
using std::string;

namespace userio {

    void aski(const char* prompt, int* iinput);
    void askr(const char* prompt, double* rinput);
    void askl(const char* prompt, bool* linput);
    void asks(const char* prompt, string* input);
    void askc(const char* prompt, char* command, string* cargs);

    void lc2uc(string* input);

    void readi(int n, int* ivar, bool* error);
    void readr(int n, double* var, bool* error);

    void getint(string input, int* a, int* n, bool* error);
    void getflt(string input, double* a, int* n, bool* error);
}

#endif //XROTOR_NOGRAPHICS_CPP_USERIO_H
