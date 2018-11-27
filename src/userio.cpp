#include <utility>

//
// Created by daniel.devries on 11/27/2018.
//
#include <algorithm>    // std::copy
#include <string>       // std::getline, std::string, std::stoi, std::stod
#include <iostream>     // std::cout, std::cin
#include <stdexcept>    // std::invalid_argument
#include "userio.h"

using namespace std;

namespace userio {

    /**
     * Presents the user with a prompt.
     *
     * @param prompt        prompt to present to the user
     * @param identifier    identified with which the actual prompt displayed will end
     */
    void _ask_helper(const char* prompt, const char* identifier) {
        cout << prompt << "    " << identifier << "> ";
    }

    /**
     * Present a prompt to the user and ask for input.
     *
     * @tparam T            known cases: int, double, string
     * @param prompt        prompt to present to the user
     * @param input         pointer to the object in which to store the result
     * @param identifier    identifier of the prompt, @see[_ask_helper(const char*, const char*)]
     */
    template <typename T>
    void _ask_helper(const char* prompt, T* input, const char* identifier) {
        _ask_helper(prompt, identifier);
        cin >> *input;
    }

    /**
     * Present a prompt to the user and ask for input with a default value.
     *
     * @tparam T            known cses: int, double, string
     * @param prompt        prompt to present to the user
     * @param input         pointer to the object in which to store the result
     * @param default_value default value in case the user does not enter anything
     * @param identifier    identifier of the prompt, @see[_ask_helper(const char*, const char*)]
     */
    template <typename T>
    void _ask_helper(const char* prompt, T* input, T default_value, const char* identifier) {
        if (*input != default_value) {
            cout << "Current value <ret takes default>: " << *input << endl;
        }
        _ask_helper(prompt, input, identifier);
    }

    /**
     * Ask the user for an integer.
     *
     * @param prompt    prompt to present to the user
     * @param iinput    pointer to the integer in which to store the result
     *
     * @note If an empty input is given, iinput will be unchanged.
     */
    void aski(const char* prompt, int* iinput) {
        _ask_helper(prompt, iinput, 999, "i");
    }

    /**
     * Ask the user for a double.
     *
     * @param prompt    prompt to present to the user
     * @param rinput    pointer to the double in which to store the result
     *
     * @note If an empty input is given, rrinput will be unchanged.
     */
    void askr(const char* prompt, double* rinput) {
        _ask_helper(prompt, rinput, 999., "r");
    }

    /**
     * Ask the user for a bool.
     *
     * @param prompt    prompt to present to the user.
     * @param linput    pointer to the bool in which to store the result
     *
     * @note Prompt will be repeated until a valid input is given.
     */
    void askl(const char* prompt, bool* linput) {
        char response = 0;
        while (response != 'Y' and response != 'N') {
            string str;
            _ask_helper(prompt, &str, "y/n");
            response = (char)toupper(str[0]);
        }
        *linput = (response == 'Y');
    }

    /**
     * Ask the user for a variable length string of characters.
     *
     * @param prompt    prompt to present to the user.
     * @param input     pointer to the string in which to store the result
     */
    void asks(const char* prompt, string* input) {
        _ask_helper(prompt, "s");
        getline(cin, *input);
    }

    /**
     * Ask the user for a command (maximum of 4 characters) and possible arguments.
     *
     * @param prompt    prompt to present to the user
     * @param command   char array with at least 4 positions
     * @param cargs     pointer to a string in which to store excess input
     *
     * @note command will be converted to upper-case.
     */
    void askc(const char* prompt, char* command, string* cargs) {
        _ask_helper(prompt, "c");

        string line;
        getline(cin, line);
        if (line.length() == 0) return;

        line.erase(0, line.find_first_not_of(" \n"));

        unsigned long k = min(line.find_first_of(" +-.,0123456789"), line.length());
        for (int i = 0; i < k; i++) {
            command[i] = (char)toupper(line[i]);
        }

        if (k < 1) k = 5;
        line = line.substr(k);
        line.erase(0, line.find_first_not_of(' '));
        line.erase(line.find_last_not_of(' ')+1, line.length());
        *cargs = line;
    }

    /**
     * Convert a string to upper-case.
     *
     * @param input     pointer to the string to convert.
     */
    void lc2uc(string* input) {
        for (char &i : *input) {
            i = (char)toupper(i);
        }
    }

    /**
     * Read a number of objects of a given type.
     *
     * @tparam T        type of the objects to read
     * @param n         number of objects to read/which are read
     * @param var       pointer to an array in which to store the objects
     * @param error     whether something went wrong
     * @param fptr      function pointer of a function which converts a string into an object of type T
     */
    template <typename T>
    void _read_helper(int n, T* var, bool* error, void (*fptr)(string, T*, int*, bool*)) {
        string line;
        getline(cin, line);

        T tmp[n];
        copy(var, var + n, tmp);

        int n_tmp = n;
        fptr(line, tmp, &n_tmp, error);
        copy(tmp, tmp + n, var);
    }

    /**
     * Read n integer values, leaving unchanged if only <return> is entered.
     *
     * @param n         number of values
     * @param ivar      array of integers
     * @param error     indicates if something went wrong while trying to read the input
     */
    void readi(int n, int* ivar, bool* error) {
        _read_helper(n, ivar, error, &getint);
    }

    /**
     * Read n double values, leaving unchanged if only <return> is entered.
     *
     * @param n         number of values
     * @param var       array of doubles
     * @param error     indicates if something went wrong while trying to read the input
     */
    void readr(int n, double* var, bool* error) {
        _read_helper(n, var, error, &getflt);
    }

    /**
     * Split string into an array of objects of a given type using given parsing function.
     *
     * @tparam T        object of which an array should be created
     * @tparam F        type of the function pointer of the parsing function
     * @param input     input string
     * @param a         resulting array of T objects
     * @param n         number of objects to expect/parsed
     * @param error     whether something went wrong
     * @param fun       function pointer of the parsing function
     *
     * @note The signature of the parsing function should be T (*)(string). This
     *       is not explicitly checked here, so this function will crash if used
     *       improperly.
     */
    template <typename T, typename F>
    void _get_helper(string input, T* a, int* n, bool* error, F fun) {
        int i = 0;
        while (input.length() > 0 && (*n != 0 && i < *n)) {
            unsigned long k = min(input.find_first_of(" ,"), input.length());
            T val;
            try {
                val = fun(input.substr(0, k));
            } catch (const invalid_argument& e) {
                *error = true;
                return;
            }
            *(a++) = val;
            input.erase(0, k + 1);
            i++;
        }
        *n = i + 1;
        *error = false;
    }

    /**
     * Parse string into an array of integers.
     *
     * Will attempt to extract no more than n numbers,
     * unless n = 0, in which case all numbers present
     * in input will be extracted.
     *
     * The resulting array will be stored in a, the
     * actual number of values extracted in n.
     *
     * @param input     input string
     * @param a         array of integers
     * @param n         number of integers
     * @param error     indicates whether something went wrong
     */
    void getint(string input, int* a, int* n, bool* error) {
        _get_helper(std::move(input), a, n, error, [](string s) -> int { return stoi(s); });
    }

    /**
     * Parse string into an array of doubles.
     *
     * Will attempt to extract no more than n numbers,
     * unless n = 0, in which case all numbers present
     * in input will be extracted.
     *
     * The resulting array will be stored in a, the
     * actual number of values extracted in n.
     *
     * @param input     input string
     * @param a         array of doubles
     * @param n         number of doubles
     * @param error     indicates whether something went wrong
     */
    void getflt(string input, double* a, int* n, bool* error) {
        _get_helper(std::move(input), a, n, error, [](string s) -> double { return stod(s); });
    }
}