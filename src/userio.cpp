//
// Created by daniel.devries on 11/27/2018.
//
#include <iostream>     // std::cout, std::cin
#include <stdexcept>    // std::invalid_argument
#include "userio.h"

using namespace std;

namespace userio {

    /**
     * Presents the user with a prompt.
     *
     * @param prompt        prompt to present to the user
     * @param identifier    string identifying what type of input is requested of the user
     */
    void _ask_helper(const string &prompt, const string &identifier) {
        cout << prompt << "    " << identifier << "> ";
    }

    /**
     * Present a prompt to the user and ask for input.
     *
     * @tparam T            known cases: int, double, string
     * @param prompt        prompt to present to the user
     * @param identifier    identifier of the prompt, @see[_ask_helper(const string &, const string &)]
     * @param input         value returned by the user
     */
    template <typename T>
    void _ask_helper(const string &prompt, const string &identifier, T &input) {
        _ask_helper(prompt, identifier);
        cin >> input;
    }

    /**
     * Present a prompt to the user and ask for input with a default value.
     *
     * @tparam T            known cases: int, double, string
     * @param prompt        prompt to present to the user
     * @param identifier    identifier of the prompt, @see[_ask_helper(const string &, const string &)]
     * @param input         value returned by the user
     * @param default_value default value in case the user does not enter anything
     */
    template <typename T>
    void _ask_helper(const string &prompt, const string &identifier, T &input, const T &default_value) {
        if (input != default_value) {
            cout << "Current value <ret takes default>: " << input << endl;
        }
        _ask_helper(prompt, identifier, input);
    }

    /**
     * Ask the user for an integer.
     *
     * @param prompt    prompt to present to the user
     * @param iinput    integer returned by the user
     *
     * @note If an empty input is given, iinput will be unchanged.
     */
    void aski(const string &prompt, int &iinput) {
        _ask_helper(prompt, "i", iinput, 999);
    }

    /**
     * Ask the user for a double.
     *
     * @param prompt    prompt to present to the user
     * @param rinput    double returned by the user
     *
     * @note If an empty input is given, rrinput will be unchanged.
     */
    void askr(const string &prompt, double &rinput) {
        _ask_helper(prompt, "r", rinput, 999.);
    }

    /**
     * Ask the user for a bool.
     *
     * @param prompt    prompt to present to the user.
     * @param linput    bool returned by the user
     *
     * @note Prompt will be repeated until a valid input is given.
     */
    void askl(const string &prompt, bool &linput) {
        char response = 0;
        while (response != 'Y' and response != 'N') {
            string str;
            _ask_helper(prompt, "y/n",  str);
            response = (char)toupper(str[0]);
        }
        linput = (response == 'Y');
    }

    /**
     * Ask the user for a variable length string of characters.
     *
     * @param prompt    prompt to present to the user.
     * @param input     string returned by the user
     */
    void asks(const string &prompt, string &input) {
        _ask_helper(prompt, "s");
        getline(cin, input);
    }

    /**
     * Ask the user for a command (maximum of 4 characters) and possible arguments.
     *
     * @param prompt    prompt to present to the user
     * @param command   command with a maximum of 4 characters
     * @param cargs     string in which to store excess input
     *
     * @note command will be converted to upper-case.
     */
    void askc(const string &prompt, string &command, string &cargs) {
        _ask_helper(prompt, "c");

        string line;
        getline(cin, line);
        if (line.length() == 0) return;

        line.erase(0, line.find_first_not_of(" \n"));

        unsigned long k = min(line.find_first_of(" +-.,0123456789"), line.length());
        command = line.substr(min(k, 4ul));
        lc2uc(command);

        if (k < 1) k = 5;
        cargs = line.substr(k);
        cargs.erase(0, cargs.find_first_not_of(' '));
        cargs.erase(cargs.find_last_not_of(' ')+1, cargs.length());
    }

    /**
     * Convert a string to upper-case.
     *
     * @param input     string to convert
     */
    void lc2uc(string &input) {
        for (auto &c : input) {
            c = (unsigned char)toupper(c);
        }
    }

    /**
     * Read a number of objects of a given type.
     *
     * @tparam T        type of the objects to read
     * @param n         number of objects to read
     * @param var       array in which to store the objects
     * @param error     whether something went wrong
     * @param fptr      either userio::getint() or userio::getflt()
     */
    template <typename T>
    void _read_helper(unsigned n, vector<T> &var, bool &error,
                      void (*fptr)(string&, vector<T>&, unsigned long&, bool&)) {
        string line;
        getline(cin, line);

        vector<T> tmp(var);
        unsigned long n_tmp = n;
        fptr(line, tmp, n_tmp, error);

        if (error) return;
        var = tmp;
    }

    /**
     * Read n integer values, leaving unchanged if only <return> is entered.
     *
     * @param n         number of values
     * @param ivar      array of integers
     * @param error     indicates if something went wrong while trying to read the input
     */
    void readi(unsigned n, ivec &ivar, bool &error) {
        _read_helper<int>(n, ivar, error, &getint);
    }

    /**
     * Read n double values, leaving unchanged if only <return> is entered.
     *
     * @param n         number of values
     * @param var       array of doubles
     * @param error     indicates if something went wrong while trying to read the input
     */
    void readr(unsigned n, vec &var, bool &error) {
        _read_helper<double>(n, var, error, &getflt);
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
    void _get_helper(string &input, vector<T> &a, unsigned long &n, bool &error, F fun) {
        a.clear();
        unsigned i = 0;
        while (input.length() > 0 && (n != 0 && i < n)) {
            unsigned long k = min(input.find_first_of(" ,"), input.length());
            try {
                T val = fun(input.substr(0, k));
                a.push_back(val);
            } catch (const invalid_argument& e) {
                error = true;
                return;
            }
            input.erase(0, k + 1);
            i++;
        }
        n = a.size();
        error = false;
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
    void getint(string &input, ivec &a, unsigned long &n, bool &error) {
        _get_helper(input, a, n, error, [](string s) -> int { return stoi(s); });
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
    void getflt(string &input, vec &a, unsigned long &n, bool &error) {
        _get_helper(input, a, n, error, [](string s) -> double { return stod(s); });
    }
}
