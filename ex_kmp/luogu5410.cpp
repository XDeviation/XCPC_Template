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
// 给两个字符串a,b，求b的exkmp和b与a的每个后缀的LCP
int main() {
    string a, b;
    cin >> a >> b;
    string c = b + '?' + a;
    vector<int> b_exkmp = exkmp(b);
    vector<int> c_exkmp = exkmp(c);
    int lenb = b.length(), lena = a.length();
    long long ans1 = lenb + 1, ans2 = c_exkmp[lenb + 1] + 1;
    for (int i = 1; i < lenb; i++) {
        ans1 ^= (i + 1LL) * (b_exkmp[i] + 1LL);
    }
    for (int i = 2; i <= lena; i++) {
        ans2 ^= i * (c_exkmp[i + lenb] + 1LL);
    }
    cout << ans1 << '\n' << ans2;
}