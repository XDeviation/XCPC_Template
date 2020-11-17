#include <bits/stdc++.h>
using namespace std;
const int maxn = 1e6 + 6;

namespace AC {
const int len = 26; //字符集个数
int tr[maxn][len], tot;
int e[maxn], fail[maxn]; // e 代表是否是结尾
void insert(string s) {
    int u = 0;
    for (int i = 0; i < s.length(); i++) {
        if (!tr[u][s[i] - 'a']) tr[u][s[i] - 'a'] = ++tot;
        u = tr[u][s[i] - 'a'];
    }
    e[u]++;
}
queue<int> q;
void build() {
    for (int i = 0; i < len; i++)
        if (tr[0][i]) q.push(tr[0][i]);
    while (q.size()) {
        int u = q.front();
        q.pop();
        for (int i = 0; i < len; i++) {
            if (tr[u][i])
                fail[tr[u][i]] = tr[fail[u]][i], q.push(tr[u][i]);
            else
                tr[u][i] = tr[fail[u]][i];
        }
    }
}
int query(string t) {
    int u = 0, res = 0;
    for (int i = 0; i < t.length(); i++) {
        u = tr[u][t[i] - 'a']; // 转移
        for (int j = u; j && e[j] != -1; j = fail[j]) {
            res += e[j], e[j] = -1;
        }
    }
    return res;
}
} // namespace AC

// 给maxn个模式串和文本串，问有多少模式串在文本串中出现

int main() {
    int n;
    cin >> n;
    string s;
    for (int i = 1; i <= n; i++)
        cin >> s, AC::insert(s);
    cin >> s;
    AC::build();
    cout << AC::query(s);
    return 0;
}