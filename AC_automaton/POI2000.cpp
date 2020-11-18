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
        if (!tr[u][s[i] - '0']) tr[u][s[i] - '0'] = ++tot;

        oi[tot] = s[i];
        u = tr[u][s[i] - '0'];
    }
    // dbg(u);
    e[u]++;
    // dbg(e);
}
void build() {
    queue<int> q;
    for (int i = 0; i < len; i++)
        if (tr[0][i]) q.push(tr[0][i]);
    while (q.size()) {
        int u = q.front();
        q.pop();
        for (int i = 0; i < len; i++) {
            if (tr[u][i])
                fail[tr[u][i]] = tr[fail[u]][i],
                e[tr[u][i]] |= e[fail[tr[u][i]]], e[tr[u][i]] |= e[u],
                q.push(tr[u][i]);
            else
                tr[u][i] = tr[fail[u]][i];
            // dbg(e);
        }
    }
    if (!tr[0][0] || !tr[0][1]) {
        puts("TAK");
        exit(0);
    }
    // dbg(e);
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
            // dbg(now, deg, i, tr[now][i]);
            dfs(tr[now][i]);
        }
    ins[now] = 0;
}
void bfs() {
    queue<int> tp;
    for (int i = 0; i <= tot; i++)
        if (!deg[i] && !e[i]) tp.push(i);
    while (!tp.empty()) {
        // dbg(tp);
        int now = tp.front();
        tp.pop();
        for (int i = 0; i < len; i++) {
            int nxtp = tr[now][i];
            if (nxtp && !e[nxtp]) {
                deg[nxtp]--;
                if (!deg[nxtp]) {
                    tp.push(nxtp);
                    dp[nxtp] = dp[now] + 1;
                    // dbg(nxtp, now);
                    nxt[nxtp] = now;
                }
            }
        }
    }
    // dbg(dp);
}
string getans() {
    vector<int> ansli;
    ansli.push_back(1);
    for (int i = 2; i <= tot; i++) {
        // dbg(i);
        if (dp[i] > dp[ansli[0]]) {
            ansli.clear();
            ansli.push_back(i);
        } else if (dp[i] == dp[ansli[0]])
            ansli.push_back(i);
    }
    // dbg(ansli);
    if (dp[ansli[0]] == 0) return "";
    vector<string> ans;
    for (auto i : ansli) {
        stack<int> tmp;
        int x = i;
        tmp.push(x);
        while (nxt[x]) {
            x = nxt[x];
            tmp.push(x);
        }
        string tmpans = "";
        while (!tmp.empty()) {
            // dbg(tmp.top());
            tmpans += oi[tmp.top()];
            tmp.pop();
        }
        ans.push_back(tmpans);
    }
    // dbg(ans);
    sort(ans.begin(), ans.end());
    return ans[0];
}
int query(string t) {
    int u = 0, res = 0;
    for (int i = 0; i < t.length(); i++) {
        u = tr[u][t[i] - '0']; // 转移
        for (int j = u; j && e[j] != -1; j = fail[j]) {
            res += e[j], e[j] = -1;
        }
    }
    return res;
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
    // for (int i = 1; i <= n; i++)
    //     cin >> inp[i];
    // sort(inp + 1, inp + 1 + n, [](string x, string y) {
    //     if (x.length() != y.length())
    //         return x.length() < y.length();
    //     else
    //         return x < y;
    // });
    // dbg(inp);
    // for (int i = 1; i <= n; i++)
    //     AC::insert(inp[i]);
    // dbg("build");
    AC::build();

    // dbg("dfs");
    AC::dfs(0);
    puts("NIE");
    return 0;
    // dbg("bfs");
    AC::bfs();
    // dbg("getans");
    cout << AC::getans() << '\n';
}