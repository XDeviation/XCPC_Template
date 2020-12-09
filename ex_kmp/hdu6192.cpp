#include <bits/stdc++.h>
using namespace std;
static auto fast_io = []() {
    std::ios::sync_with_stdio(false); // turn off sync
    cin.tie(nullptr);                 // untie in/out streams
    return 0;
}();
vector<int> exkmp(string s) {
    int n = (int) s.length();
    vector<int> z(n);
    for (int i = 1, l = 0, r = 0; i < n; ++i) {
        if (i <= r) z[i] = min(r - i + 1, z[i - l]);
        while (i + z[i] < n && s[z[i]] == s[i + z[i]])
            ++z[i];
        if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
    }
    return z;
}
const int maxn = 1e5 + 7;
const unsigned int mi = 131;
unsigned int ex[maxn], has[maxn];
string s;
unsigned int _seg(int l, int r, int k, char x) {
    return has[r] - (l - 1 >= 0 ? has[l - 1] : 0) * ex[r - l + 1] -
           (k >= l && k <= r ? ex[r - k] * s[k] - ex[r - k] * x : 0);
};
// 修改字符串中一个字符，使得循环节最短（最后一个可以不取完全）
int main() {
    int n;
    while (cin >> n >> s) {
        // dbg(n, s);
        vector<int> s_exkmp = exkmp(s);
        ex[0] = 1, has[0] = s[0];
        // dbg(s_exkmp);
        for (int i = 1; i <= n; i++)
            ex[i] = ex[i - 1] * mi;
        for (int i = 1; i < n; i++)
            has[i] = has[i - 1] * mi + s[i];
        for (int i = 1; i < n; i++) {
            int v = s_exkmp[i];
            if (v + i == n) {
                cout << i << ' ' << n << '\n';
                break;
            } else {
                int p = 0;
                // dbg(i);
                p += (_seg(0, n - 1 - i, v, s[i + v]) ==
                      _seg(i, n - 1, v, s[i + v]));
                p += (_seg(0, n - 1 - i, i + v, s[v]) ==
                      _seg(i, n - 1, i + v, s[v]));
                if (p) {
                    cout << i << ' ' << p << '\n';
                    break;
                }
            }
        }
        if (n == 1) cout << "1 1\n";
    }
}