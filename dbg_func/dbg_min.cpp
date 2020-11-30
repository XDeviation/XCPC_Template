#include <bits/stdc++.h>
using namespace std;

#define show(arg) cerr << "\033[36m" << #arg << "\033[0m = " << arg << endl;
#define dbg(...) debug(string(#__VA_ARGS__).begin(), __VA_ARGS__)

void debug(string::iterator it) {
    cerr << endl;
}

void debug(string::iterator it, auto a, auto... args) {
    cerr << "\033[36m";
    int cnt = 0;
    while (*it != ',' || cnt) {
        if (*it == '(' || *it == '{') cnt++;
        if (*it == ')' || *it == '}') cnt--;
        if (*it != ' ') cerr << *it;
        ++it;
    }
    cerr << "\033[0m"
         << " = " << a << "   ";
    debug(++it, args...);
}