#include <bits/stdc++.h>
// #include <dbg_func>
using namespace std;
const int maxn = 1e6 + 6;
vector<pair<int, int>> G[maxn], R[maxn];
namespace AC {
const int len = 2; //字符集个数
int tr[maxn][len], tot;
int e[maxn], fail[maxn];
int vis[maxn], ins[maxn];
int deg[maxn], ord[maxn], dp[maxn];
int mxd[maxn];
char oi[maxn];
vector<int> pi;
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
                fail[v] = tr[fail[u]][i], e[v] |= e[fail[v]], q.push(v);
            else
                v = tr[fail[u]][i];
        }
    }
}
void dfs(int now) {
    vis[now] = 1;
    pi.push_back(now);
    for (int i = 0; i < len; i++) {
        int v = tr[now][i];
        if (!vis[v] && !e[v]) dfs(v);
    }
}

int topo() {
    tot = 0;
    for (auto u : pi) {
        for (int i = 0; i < len; i++) {
            int v = tr[u][i];
            if (vis[v]) {
                G[u].push_back({i, v}), deg[v]++;
            }
        }
    }
    queue<int> q;
    for (auto it : pi)
        if (!deg[it]) q.push(it);
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        ord[++tot] = u;
        for (auto it : G[u]) {
            if (--deg[it.second] == 0) q.push(it.second);
        }
    }
    return tot == pi.size();
}
void dfs2(int u) {
    mxd[u] = dp[u];
    for (auto it : R[u]) {
        dfs2(it.second);
        mxd[u] = max(mxd[u], mxd[it.second]);
    }
}
void print(int u) {
    for (auto it : R[u]) {
        if (mxd[it.second] == mxd[u]) {
            cout << it.first;
            print(it.second);
            break;
        }
    }
}
void solve() {
    for (int i = 1; i <= tot; i++) {
        int u = ord[i];
        for (auto it : G[u]) {
            dp[it.second] = max(dp[it.second], dp[u] + 1);
        }
    }
    int ret = 0;
    for (int i = 1; i <= tot; i++)
        ret = max(ret, dp[i]);
    for (int i = 1; i <= tot; i++) {
        int u = ord[i];
        for (auto it : G[u]) {
            if (dp[it.second] == dp[u] + 1) {
                R[u].push_back(it);
            }
        }
        sort(R[u].begin(), R[u].end());
    }
    dfs2(0);
    print(0);
}
} // namespace AC

// 2016 ICPC 香港 K 题 trie 图上找最长链
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
    if (!AC::topo()) {
        cout << "-1\n";
        return 0;
    }
    AC::solve();
    return 0;
}