/*
 * Codeforces Contest
 * Au: Lost_Deviation
 */
#include <bits/stdc++.h>
#ifndef open_dbg_func
#define dbg(args...) (args)
#endif
using namespace std;
#define ll long long
#define bll unsigned long long
const int maxn = 1e5 + 7;
static auto x = []() {
    std::ios::sync_with_stdio(false); // turn off sync
    cin.tie(nullptr);                 // untie in/out streams
    return 0;
}();
vector<string> split(string str, const string &delimiters) {
    regex del(delimiters);
    vector<string> v(sregex_token_iterator(str.begin(), str.end(), del, -1),
                     sregex_token_iterator());
    return v;
}
struct T {
    int x, y;
};
int cmp(const T &a, const T &b) {
    if (a.x != b.x)
        return (a.x < b.x);
    else
        return (a.y < b.y);
}
int main(void) {
    return 0;
}
