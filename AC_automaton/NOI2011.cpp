#include <bits/stdc++.h>
using namespace std;
const int maxn = 1e5 + 6;

namespace AC {
const int len = 26; //字符集个数
int tr[maxn][len], tot;
int fail[maxn], last[maxn], fa[maxn], siz[maxn], dfsq[maxn];
vector<int> failtree[maxn];
int enode[maxn]; // 第i个串对应的结束节点编号
vector<int> e[maxn];
vector<int> nxt[maxn];
vector<pair<int, int>> que[maxn]; // first id second x
int ans[maxn];
int n, cnt, m; // dfs序
int tre[maxn]; //树状数组

void insert(string s) {
    int u = 0;
    for (int i = 0; i < s.length(); i++) {
        if (s[i] == 'B') {
            u = fa[u];
        } else if (s[i] == 'P') {
            e[u].push_back(++n);
            enode[n] = u;
        } else {
            int &v = tr[u][s[i] - 'a'];
            if (!v) {
                v = ++tot;
                nxt[u].push_back(tot);
                fa[tot] = u;
            }
            u = v;
        }
    }
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
                fail[v] = tr[fail[u]][i], q.push(v);
            else
                v = tr[fail[u]][i];
        }
    }
    for (int i = 1; i <= tot; i++)
        failtree[fail[i]].push_back(i);
}

void dfs(int u) {
    siz[u] = 1;
    dfsq[u] = ++cnt;
    for (auto i : failtree[u]) {
        dfs(i);
        siz[u] += siz[i];
    }
    // dbg(dfsq, siz);
}

void readquestion() {
    cin >> m;
    for (int i = 1; i <= m; i++) {
        int x, y;
        cin >> x >> y;
        que[y].push_back({i, x});
    }
}

void treeinsert(int x, int inc) {
    for (x; x <= cnt; x += (x & -x))
        tre[x] += inc;
}

int sum(int x) {
    int res = 0;
    for (x; x; x -= (x & -x))
        res += tre[x];
    return res;
}

void work(int u) {
    treeinsert(dfsq[u], 1);
    for (auto it : e[u]) {
        for (auto jt : que[it]) {
            int v = enode[jt.second];
            ans[jt.first] += sum(dfsq[v] + siz[v] - 1) - sum(dfsq[v] - 1);
        }
    }
    for (auto it : nxt[u])
        work(it);
    treeinsert(dfsq[u], -1);
}

void print() {
    for (int i = 1; i <= m; i++)
        cout << ans[i] << '\n';
}
} // namespace AC

// 阿狸的打字机

int main() {
    string s;
    cin >> s;
    AC::insert(s);
    AC::build();
    AC::dfs(0);
    AC::readquestion();
    AC::work(0);
    AC::print();
    return 0;
}