#include <bits/stdc++.h>
// #include <dbg_func>
using namespace std;
const int maxn = 1e6 + 6;
namespace AC {
const int len = 2; //字符集个数
int tr[maxn][len], tot;
int e[maxn], fail[maxn];
int vis[maxn], ins[maxn];
int deg[maxn], nxt[maxn], dp[maxn];
char oi[maxn];
void insert(string s) {
    int u = 0;
    for (int i = 0; i < s.length(); i++) {
        int &v = tr[u][s[i] - '0'];
        if (!v) v = ++tot;
        oi[tot] = s[i];
        u = v;
    }
    e[u]++;
}
void build() {
    queue<int> q;
    for (int i = 0; i < len; i++)
        if (tr[0][i]) q.push(tr[0][i]);
    while (q.size()) {
        int u = q.front();
        q.pop();
        for (int i = 0; i < len; i++) {
            int &v = tr[u][i];
            if (v)
                fail[v] = tr[fail[u]][i], e[v] |= e[fail[v]], e[v] |= e[u],
                q.push(v);
            else
                v = tr[fail[u]][i];
            // dbg(e);
        }
    }
    if (!tr[0][0] || !tr[0][1]) {
        puts("TAK");
        exit(0);
    }
}
void dfs(int now) {
    if (ins[now]) {
        puts("TAK");
        exit(0);
    }
    if (vis[now] || e[now]) return;
    vis[now] = ins[now] = 1;
    for (int i = 0; i < len; i++)
        if (tr[now][i]) {
            deg[tr[now][i]]++;
            dfs(tr[now][i]);
        }
    ins[now] = 0;
}
} // namespace AC

// 给maxn个模式串和文本串，问有多少模式串在文本串中出现
string inp[maxn];
int main() {
    int n;
    cin >> n;
    string s;
    for (int i = 1; i <= n; i++) {
        cin >> s;
        AC::insert(s);
    }
    AC::build();
    AC::dfs(0);
    puts("NIE");
    return 0;
}