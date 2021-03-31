 
```cpp
//./suffix_automaton_plus/luogu6139.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	#define MAXN 2000000 // 双倍字符串长度
     4	#define CHAR_NUM 30  // 字符集个数，注意修改下方的 (-'a')
     5	struct exSAM {
     6	    int len[MAXN];            // 节点长度
     7	    int link[MAXN];           // 后缀链接，link
     8	    int next[MAXN][CHAR_NUM]; // 转移
     9	    int tot;                  // 节点总数：[0, tot)
    10	    void init() {
    11	        tot = 1;
    12	        link[0] = -1;
    13	    }
    14	    int insertSAM(int last, int c) {
    15	        int cur = next[last][c];
    16	        if (len[cur]) return cur;
    17	        len[cur] = len[last] + 1;
    18	        int p = link[last];
    19	        while (p != -1) {
    20	            if (!next[p][c])
    21	                next[p][c] = cur;
    22	            else
    23	                break;
    24	            p = link[p];
    25	        }
    26	        if (p == -1) {
    27	            link[cur] = 0;
    28	            return cur;
    29	        }
    30	        int q = next[p][c];
    31	        if (len[p] + 1 == len[q]) {
    32	            link[cur] = q;
    33	            return cur;
    34	        }
    35	        int clone = tot++;
    36	        for (int i = 0; i < CHAR_NUM; ++i)
    37	            next[clone][i] = len[next[q][i]] != 0 ? next[q][i] : 0;
    38	        len[clone] = len[p] + 1;
    39	        while (p != -1 && next[p][c] == q) {
    40	            next[p][c] = clone;
    41	            p = link[p];
    42	        }
    43	        link[clone] = link[q];
    44	        link[cur] = clone;
    45	        link[q] = clone;
    46	        return cur;
    47	    }
    48	    int insertTrie(int cur, int c) {
    49	        if (next[cur][c]) return next[cur][c];
    50	        return next[cur][c] = tot++;
    51	    }
    52	    void insert(const string &s) {
    53	        int root = 0;
    54	        for (auto ch : s)
    55	            root = insertTrie(root, ch - 'a');
    56	    }
    57	    void insert(const char *s, int n) {
    58	        int root = 0;
    59	        for (int i = 0; i < n; ++i)
    60	            root = insertTrie(root, s[i] - 'a');
    61	    }
    62	    void build() {
    63	        queue<pair<int, int>> q;
    64	        for (int i = 0; i < 26; ++i)
    65	            if (next[0][i]) q.push({i, 0});
    66	        while (!q.empty()) {
    67	            auto item = q.front();
    68	            q.pop();
    69	            auto last = insertSAM(item.second, item.first);
    70	            for (int i = 0; i < 26; ++i)
    71	                if (next[last][i]) q.push({i, last});
    72	        }
    73	    }
    74	} exSam;
    75	char s[1000100];

    76	// 给定 n 个由小写字母组成的字符串
    77	// s1,s2 … sn​，求本质不同的子串个数。（不包含空串）

    78	int main() {
    79	    int n;
    80	    cin >> n;
    81	    exSam.init();
    82	    for (int i = 0; i < n; ++i) {
    83	        cin >> s;
    84	        int len = strlen(s);
    85	        exSam.insert(s, len);
    86	    }
    87	    exSam.build();
    88	    long long ans = 0;
    89	    for (int i = 1; i < exSam.tot; ++i) {
    90	        ans += exSam.len[i] - exSam.len[exSam.link[i]];
    91	    }
    92	    cout << ans << endl;
    93	}
 
```
 
```cpp
//./suffix_automaton_plus/spoj.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	#define MAXN 2000000 // 双倍字符串长度
     4	#define CHAR_NUM 30  // 字符集个数，注意修改下方的 (-'a')
     5	#define NUM 15       // 字符串个数
     6	struct exSAM {
     7	    int len[MAXN];            // 节点长度
     8	    int link[MAXN];           // 后缀链接，link
     9	    int next[MAXN][CHAR_NUM]; // 转移
    10	    int tot;                  // 节点总数：[0, tot)
    11	    int lenSorted[MAXN];  // 按照 len 排序后的数组，仅排序 [1, tot)
    12	                          // 部分，最终下标范围 [0, tot - 1)
    13	    int sizeC[MAXN][NUM]; // 表示某个字符串的子串个数
    14	    int curString;        // 字符串实际个数
    15	    /**
    16	     * 计数排序使用的辅助空间数组
    17	     */
    18	    int lc[MAXN]; // 统计个数
    19	    void init() {
    20	        tot = 1;
    21	        link[0] = -1;
    22	    }
    23	    int insertSAM(int last, int c) {
    24	        int cur = next[last][c];
    25	        len[cur] = len[last] + 1;
    26	        int p = link[last];
    27	        while (p != -1) {
    28	            if (!next[p][c])
    29	                next[p][c] = cur;
    30	            else
    31	                break;
    32	            p = link[p];
    33	        }
    34	        if (p == -1) {
    35	            link[cur] = 0;
    36	            return cur;
    37	        }
    38	        int q = next[p][c];
    39	        if (len[p] + 1 == len[q]) {
    40	            link[cur] = q;
    41	            return cur;
    42	        }
    43	        int clone = tot++;
    44	        for (int i = 0; i < CHAR_NUM; ++i)
    45	            next[clone][i] = len[next[q][i]] != 0 ? next[q][i] : 0;
    46	        len[clone] = len[p] + 1;
    47	        while (p != -1 && next[p][c] == q) {
    48	            next[p][c] = clone;
    49	            p = link[p];
    50	        }
    51	        link[clone] = link[q];
    52	        link[cur] = clone;
    53	        link[q] = clone;
    54	        return cur;
    55	    }
    56	    int insertTrie(int cur, int c) {
    57	        if (!next[cur][c]) next[cur][c] = tot++;
    58	        sizeC[next[cur][c]][curString]++;
    59	        return next[cur][c];
    60	    }
    61	    void insert(const string &s) {
    62	        int root = 0;
    63	        for (auto ch : s)
    64	            root = insertTrie(root, ch - 'a');
    65	        curString++;
    66	    }
    67	    void insert(const char *s, int n) {
    68	        int root = 0;
    69	        for (int i = 0; i < n; ++i)
    70	            root = insertTrie(root, s[i] - 'a');
    71	        curString++;
    72	    }
    73	    void build() {
    74	        queue<pair<int, int>> q;
    75	        for (int i = 0; i < 26; ++i)
    76	            if (next[0][i]) q.push({i, 0});
    77	        while (!q.empty()) {
    78	            auto item = q.front();
    79	            q.pop();
    80	            auto last = insertSAM(item.second, item.first);
    81	            for (int i = 0; i < 26; ++i)
    82	                if (next[last][i]) q.push({i, last});
    83	        }
    84	    }
    85	    void sortLen() {
    86	        for (int i = 1; i < tot; ++i)
    87	            lc[i] = 0;
    88	        for (int i = 1; i < tot; ++i)
    89	            lc[len[i]]++;
    90	        for (int i = 2; i < tot; ++i)
    91	            lc[i] += lc[i - 1];
    92	        for (int i = 1; i < tot; ++i)
    93	            lenSorted[--lc[len[i]]] = i;
    94	    }
    95	    void getSizeLen() {
    96	        for (int i = tot - 2; i >= 0; --i)
    97	            for (int j = 0; j < curString; ++j)
    98	                sizeC[link[lenSorted[i]]][j] += sizeC[lenSorted[i]][j];
    99	    }
   100	    void debug() {
   101	        cout << "     i      len      link       ";
   102	        for (int i = 0; i < 26; ++i)
   103	            cout << "  " << (char) ('a' + i);
   104	        cout << endl;
   105	        for (int i = 0; i < tot; ++i) {
   106	            cout << "i: " << setw(3) << i << " len: " << setw(3) << len[i]
   107	                 << " link: " << setw(3) << link[i] << " Next: ";
   108	            for (int j = 0; j < CHAR_NUM; ++j) {
   109	                cout << setw(3) << next[i][j];
   110	            }
   111	            cout << endl;
   112	        }
   113	    }
   114	} exSam;
   115	// 找所有字符串的最长公共子串
   116	int main() {
   117	    exSam.init();
   118	    string s;
   119	    while (cin >> s)
   120	        exSam.insert(s);
   121	    exSam.build();
   122	    exSam.sortLen();
   123	    exSam.getSizeLen();
   124	    int ans = 0;
   125	    for (int i = 0; i < exSam.tot; ++i) {
   126	        bool flag = true;
   127	        for (int j = 0; j < exSam.curString; ++j) {
   128	            if (!exSam.sizeC[i][j]) {
   129	                flag = false;
   130	                break;
   131	            }
   132	        }
   133	        if (flag) ans = max(ans, exSam.len[i]);
   134	    }
   135	    cout << ans << endl;
   136	}

 
```
 
```cpp
//./suffix_tree/luogu3804.cpp
     1	/*

     2	   在这份代码实现中，rem表示的是len，0号节点的实际编号为1

     3	   上文中提到的未插入后缀数rem没有被维护，而是被len和now表示

     4	   len[u]和start[u]表示u的父边在s中的起点和长度，
     5	   即u的父边代表s[start[u] ... start[u]+len[u]-1]

     6	    求字符串中所有子串数量（大于2） 注意字符串从1开始！
     7	 */
     8	#include <bits/stdc++.h>
     9	using namespace std;
    10	#define ll long long
    11	const int maxn = 1e5 + 7;
    12	const int inf = 1e9;
    13	int siz[maxn];
    14	struct suffixTree {
    15	    int link[maxn], len[maxn], start[maxn], s[maxn], n, tail, now, rem;
    16	    map<int, int> ch[maxn];
    17	    suffixTree()
    18	        : tail(1)
    19	        , n(0)
    20	        , rem(0)
    21	        , now(1) {
    22	        len[0] = inf;
    23	    }
    24	    int newnode(int st, int le) {
    25	        link[++tail] = 1;
    26	        start[tail] = st;
    27	        len[tail] = le;
    28	        return tail;
    29	    }
    30	    void extend(int x) {
    31	        s[++n] = x;
    32	        rem++;
    33	        for (int last = 1; rem;) {
    34	            while (rem > len[ch[now][s[n - rem + 1]]]) {
    35	                rem -= len[now = ch[now][s[n - rem + 1]]];
    36	            }
    37	            int &v = ch[now][s[n - rem + 1]];
    38	            int c = s[start[v] + rem - 1];
    39	            if (!v || x == c) {
    40	                link[last] = now;
    41	                last = now;
    42	                if (!v)
    43	                    v = newnode(n, inf);
    44	                else
    45	                    break;
    46	            } else {
    47	                int u = newnode(start[v], rem - 1);
    48	                ch[u][c] = v;
    49	                ch[u][x] = newnode(n, inf);
    50	                start[v] += rem - 1;
    51	                len[v] -= rem - 1;
    52	                link[last] = v = u;
    53	                last = u;
    54	            }
    55	            if (now == 1)
    56	                rem--;
    57	            else
    58	                now = link[now];
    59	        }
    60	    }
    61	} sft;
    62	ll ans = 0;
    63	int dfs(int u, int depth) {
    64	    if (depth >= inf) return 1;
    65	    siz[u] = 0;
    66	    for (auto it : sft.ch[u]) {
    67	        int d = dfs(it.second, depth + sft.len[it.second]);
    68	        siz[u] += d;
    69	    }
    70	    if (siz[u] >= 2) ans = max(ans, (ll) siz[u] * (ll) depth);
    71	    return siz[u];
    72	}
    73	char s[maxn];
    74	int main() {
    75	    scanf("%s", s + 1);
    76	    int n = std::strlen(s + 1);

    77	    for (int i = 1; i <= n; ++i)
    78	        sft.extend(s[i] - 'a');
    79	    sft.extend(26);
    80	    dfs(1, 0);
    81	    printf("%lld", ans);
    82	    return 0;
    83	}
 
```
 
```cpp
//./3D_convex/luogu4724.cpp
     1	#include <cmath>
     2	#include <cstdio>
     3	#include <cstdlib>
     4	#include <iostream>
     5	using namespace std;
     6	const int N = 2010;
     7	const double eps = 1e-9;
     8	int n, cnt, vis[N][N];
     9	double ans;
    10	double Rand() {
    11	    return rand() / (double) RAND_MAX;
    12	}
    13	double reps() {
    14	    return (Rand() - 0.5) * eps;
    15	}
    16	struct Node {
    17	    double x, y, z;
    18	    void shake() {
    19	        x += reps();
    20	        y += reps();
    21	        z += reps();
    22	    }
    23	    double len() {
    24	        return sqrt(x * x + y * y + z * z);
    25	    }
    26	    Node operator-(Node A) {
    27	        return (Node){x - A.x, y - A.y, z - A.z};
    28	    }
    29	    Node operator*(Node A) {
    30	        return (Node){y * A.z - z * A.y, z * A.x - x * A.z, x * A.y - y * A.x};
    31	    }
    32	    double operator&(Node A) {
    33	        return x * A.x + y * A.y + z * A.z;
    34	    }
    35	} A[N];
    36	struct Face {
    37	    int v[3];
    38	    Node Normal() {
    39	        return (A[v[1]] - A[v[0]]) * (A[v[2]] - A[v[0]]);
    40	    }
    41	    double area() {
    42	        return Normal().len() / 2.0;
    43	    }
    44	} f[N], C[N];
    45	int see(Face a, Node b) {
    46	    return ((b - A[a.v[0]]) & a.Normal()) > 0;
    47	}
    48	void Convex_3D() {
    49	    f[++cnt] = (Face){1, 2, 3};
    50	    f[++cnt] = (Face){3, 2, 1};
    51	    for (int i = 4, cc = 0; i <= n; i++) {
    52	        for (int j = 1, v; j <= cnt; j++) {
    53	            if (!(v = see(f[j], A[i]))) C[++cc] = f[j];
    54	            for (int k = 0; k < 3; k++)
    55	                vis[f[j].v[k]][f[j].v[(k + 1) % 3]] = v;
    56	        }
    57	        for (int j = 1; j <= cnt; j++)
    58	            for (int k = 0; k < 3; k++) {
    59	                int x = f[j].v[k], y = f[j].v[(k + 1) % 3];
    60	                if (vis[x][y] && !vis[y][x]) C[++cc] = (Face){x, y, i};
    61	            }
    62	        for (int j = 1; j <= cc; j++)
    63	            f[j] = C[j];
    64	        cnt = cc;
    65	        cc = 0;
    66	    }
    67	}
    68	// 给空间n个点 求凸包表面积
    69	int main() {
    70	    cin >> n;
    71	    for (int i = 1; i <= n; i++)
    72	        cin >> A[i].x >> A[i].y >> A[i].z, A[i].shake();
    73	    Convex_3D();
    74	    for (int i = 1; i <= cnt; i++)
    75	        ans += f[i].area();
    76	    printf("%.3f\n", ans);
    77	}
 
```
 
```cpp
//./LCS/lcs.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;

     3	typedef long long ll;
     4	const int BIT = 808, E = 62;
     5	int n1, n2, m, h, i, j, p, ans;
     6	int s[50000], t[50000];
     7	struct Num {
     8	    ll x[BIT];
     9	    Num() {
    10	        for (int i = 0; i < BIT; i++)
    11	            x[i] = 0;
    12	    }
    13	    void set(int p) {
    14	        x[p / E] |= 1LL << (p % E);
    15	    }
    16	    Num operator&(Num b) {
    17	        Num c;
    18	        for (int i = 0; i <= m; i++)
    19	            c.x[i] = x[i] & b.x[i];
    20	        return c;
    21	    }
    22	    Num operator|(Num b) {
    23	        Num c;
    24	        for (int i = 0; i <= m; i++)
    25	            c.x[i] = x[i] | b.x[i];
    26	        return c;
    27	    }
    28	    Num operator^(Num b) {
    29	        Num c;
    30	        for (int i = 0; i <= m; i++)
    31	            c.x[i] = x[i] ^ b.x[i];
    32	        return c;
    33	    }
    34	    Num operator-(Num b) {
    35	        Num c;
    36	        for (int i = 0; i <= m; i++)
    37	            c.x[i] = x[i] - b.x[i];
    38	        for (int i = 0; i < m; i++)
    39	            if (c.x[i] < 0) c.x[i] += (1LL << E), c.x[i + 1];
    40	        return c;
    41	    }
    42	    void shl() {
    43	        ll o = 1, p;
    44	        for (int i = 0; i <= m; o = p, i++) {
    45	            p = x[i] & (1LL << h), (x[i] <<= 1) &= ~(1LL << (h + 1));
    46	            if (o) x[i] |= 1;
    47	        }
    48	    }
    49	} ap[100007], x, row[2];
    50	int my_hash(int x) {
    51	    if (x == 'A') return 0;
    52	    if (x == 'C') return 1;
    53	    if (x == 'G') return 2;
    54	    return 3;
    55	}
    56	int main() {
    57	    // scanf("%d%d%s%s", &n1, &n2, s, t);
    58	    int n;
    59	    cin >> n;
    60	    n1 = n, n2 = n;
    61	    for (int it = 1; it <= n; it++)
    62	        cin >> s[it];
    63	    for (int it = 1; it <= n; it++)
    64	        cin >> t[it];
    65	    i = 0;
    66	    for (m = (n1 - 1) / E, h = (m ? E : n1) - 1; i < n1; i++)
    67	        ap[s[i]].set(i);
    68	    for (i = 0; i < n2; i++) {
    69	        p ^= 1;
    70	        x = ap[t[i]] | row[p ^ 1];
    71	        row[p ^ 1].shl();
    72	        row[p] = x & ((x - row[p ^ 1]) ^ x);
    73	    }
    74	    for (i = m; ~i; i--)
    75	        for (j = h; ~j; j--)
    76	            if (row[p].x[i] & (1LL << j)) ans++;
    77	    printf("%d", ans);
    78	}
 
```
 
```cpp
//./pb_ds_queue/luogu1552.cpp
     1	#include <bits/extc++.h>
     2	#include <bits/stdc++.h>
     3	#ifndef ONLINE_JUDGE
     4	#include "dbg_func"
     5	#else
     6	#define dbg(...) (__VA_ARGS__)
     7	#endif
     8	/*
     9	若编译器没有 extc++.h 使用如下三个库
    10	#include <ext/pb_ds/assoc_container.hpp>
    11	#include <ext/pb_ds/priority_queue.hpp>
    12	#include <ext/pb_ds/tree_policy.hpp>
    13	*/
    14	#define param int, less<int>, pairing_heap_tag
    15	using namespace std;
    16	using namespace __gnu_pbds;
    17	const int maxn = 1e5 + 7;
    18	#define int long long
    19	/*
    20	__gnu_pbds::priority_queue<int, greater<int>> que1;                       // （默认为配对堆）
    21	__gnu_pbds::priority_queue<int, less<int>> que1;                          // （默认为配对堆）
    22	__gnu_pbds::priority_queue<int, greater<int>, pairing_heap_tag> que2;     // 配对堆（应该是最快的）
    23	__gnu_pbds::priority_queue<int, greater<int>, binary_heap_tag> que3;      // 二叉堆
    24	__gnu_pbds::priority_queue<int, greater<int>, binomial_heap_tag> que4;    // 二项堆
    25	__gnu_pbds::priority_queue<int, greater<int>, rc_binomial_heap_tag> que5; //
    26	__gnu_pbds::priority_queue<int, greater<int>, thin_heap_tag> que6;        // 斐波那契堆
    27	只有 push pop join 用 binary_heap
    28	有 modify 用 pairing_heap & thin_heap <- thin_heap 可能会MLE，慎用
    29	优化 Dijkstra 用 pairing_heap
    30	std 垃圾
    31	less 大根堆 greater 小根堆
    32	例题 猴子选国王 并查集 + 左偏树
    33	*/
    34	__gnu_pbds::priority_queue<param> pq[maxn];

    35	struct DisjointSet {
    36	    int p[maxn];
    37	    void init(int n) {
    38	        for (int i = 1; i <= n; i++) {
    39	            p[i] = i;
    40	        }
    41	    }
    42	    int find(int u) { return u == p[u] ? u : p[u] = find(p[u]); }
    43	    bool isJoined(int u, int v) { return find(u) == find(v); }

    44	    void Union(int u, int v) {
    45	        u = find(u);
    46	        v = find(v);
    47	        p[v] = u;
    48	    }
    49	} s;
    50	int c[maxn], l[maxn], deg[maxn], totc[maxn];

    51	signed main() {
    52	    /*
    53	    修改和删除示例
    54	    __gnu_pbds::priority_queue<int> testplus;
    55	    auto it = testplus.push(0);
    56	    testplus.push(1);
    57	    testplus.push(2);
    58	    testplus.modify(it, 3);
    59	    cout << testplus.top();
    60	    testplus.erase(it);
    61	    cout << testplus.top();
    62	    */
    63	    int n, m;
    64	    cin >> n >> m;
    65	    s.init(n);
    66	    for (int i = 1; i <= n; i++) {
    67	        int bi;
    68	        cin >> bi >> c[i] >> l[i];
    69	        s.p[i] = bi;
    70	        deg[bi]++;
    71	        totc[i] = c[i];
    72	        pq[i].push(c[i]);
    73	    }
    74	    queue<int> ninjya;
    75	    for (int i = 1; i <= n; i++)
    76	        if (deg[i] == 0) {
    77	            ninjya.push(i);
    78	        }
    79	    int ans = 0;
    80	    while (ninjya.size()) {
    81	        int now = ninjya.front();
    82	        ninjya.pop();
    83	        while (totc[now] > m) {
    84	            totc[now] -= pq[now].top();
    85	            pq[now].pop();
    86	        }
    87	        // cout << now << ' ' << totc[now] << endl;
    88	        dbg(now, totc[now], l[now], pq[now].size());
    89	        ans = max(ans, l[now] * (long long)pq[now].size());
    90	        pq[s.p[now]].join(pq[now]);
    91	        totc[s.p[now]] += totc[now];
    92	        deg[s.p[now]]--;
    93	        if (deg[s.p[now]] == 0)
    94	            ninjya.push(s.p[now]);
    95	    }
    96	    cout << ans << endl;
    97	    return 0;
    98	}
 
```
 
```cpp
//./pb_ds_queue/pb_ds_queue.cpp
     1	#include <bits/extc++.h>
     2	#include <bits/stdc++.h>
     3	/*
     4	若编译器没有 extc++.h 使用如下三个库
     5	#include <ext/pb_ds/assoc_container.hpp>
     6	#include <ext/pb_ds/priority_queue.hpp>
     7	#include <ext/pb_ds/tree_policy.hpp>
     8	*/
     9	#define param int, less<int>, pairing_heap_tag
    10	using namespace std;
    11	using namespace __gnu_pbds;
    12	const int maxn = 1e5 + 7;
    13	/*
    14	__gnu_pbds::priority_queue<int, greater<int>> que1;                       // （默认为配对堆）
    15	__gnu_pbds::priority_queue<int, less<int>> que1;                          // （默认为配对堆）
    16	__gnu_pbds::priority_queue<int, greater<int>, pairing_heap_tag> que2;     // 配对堆（应该是最快的）
    17	__gnu_pbds::priority_queue<int, greater<int>, binary_heap_tag> que3;      // 二叉堆
    18	__gnu_pbds::priority_queue<int, greater<int>, binomial_heap_tag> que4;    // 二项堆
    19	__gnu_pbds::priority_queue<int, greater<int>, rc_binomial_heap_tag> que5; //
    20	__gnu_pbds::priority_queue<int, greater<int>, thin_heap_tag> que6;        // 斐波那契堆
    21	只有 push pop join 用 binary_heap
    22	有 modify 用 pairing_heap & thin_heap <- thin_heap 可能会MLE，慎用
    23	优化 Dijkstra 用 pairing_heap
    24	std 垃圾
    25	less 大根堆 greater 小根堆
    26	例题 猴子选国王 并查集 + 左偏树
    27	*/
    28	__gnu_pbds::priority_queue<param> qs[maxn];

    29	struct DisjointSet {
    30	    int p[maxn];
    31	    void init(int n) {
    32	        for (int i = 1; i <= n; i++) {
    33	            p[i] = i;
    34	        }
    35	    }
    36	    int find(int u) { return u == p[u] ? u : p[u] = find(p[u]); }
    37	    bool isJoined(int u, int v) { return find(u) == find(v); }

    38	    void Union(int u, int v) {
    39	        u = find(u);
    40	        v = find(v);
    41	        p[v] = u;
    42	    }
    43	} s;

    44	void jointHeap(int u, int v) {
    45	    u = s.find(u);
    46	    v = s.find(v);
    47	    int utop = qs[u].top();
    48	    qs[u].pop();
    49	    qs[u].push(utop / 2);
    50	    int vtop = qs[v].top();
    51	    qs[v].pop();
    52	    qs[v].push(vtop / 2);
    53	    qs[u].join(qs[v]);
    54	}

    55	int main() {
    56	    /*
    57	    修改和删除示例
    58	    __gnu_pbds::priority_queue<int> testplus;
    59	    auto it = testplus.push(0);
    60	    testplus.push(1);
    61	    testplus.push(2);
    62	    testplus.modify(it, 3);
    63	    cout << testplus.top();
    64	    testplus.erase(it);
    65	    cout << testplus.top();
    66	    */
    67	    int n, m;
    68	    while (cin >> n) {
    69	        s.init(n);
    70	        for (int i = 1; i <= n; i++) {
    71	            int val;
    72	            qs[i].clear();
    73	            cin >> val;
    74	            qs[i].push(val);
    75	        }
    76	        cin >> m;
    77	        while (m--) {
    78	            int u, v;
    79	            cin >> u >> v;
    80	            if (s.isJoined(u, v)) {
    81	                cout << -1 << endl;
    82	                continue;
    83	            }
    84	            jointHeap(u, v);
    85	            cout << qs[s.find(u)].top() << endl;
    86	            s.Union(u, v);
    87	        }
    88	    }
    89	    return 0;
    90	}
 
```
 
```cpp
//./pb_ds_queue/luogu3377.cpp
     1	#include <bits/extc++.h>
     2	#include <bits/stdc++.h>
     3	/*
     4	若编译器没有 extc++.h 使用如下三个库
     5	#include <ext/pb_ds/assoc_container.hpp>
     6	#include <ext/pb_ds/priority_queue.hpp>
     7	#include <ext/pb_ds/tree_policy.hpp>
     8	*/
     9	class T {

    10	  public:
    11	    int val, rnk;
    12	    T(int val, int rnk) {
    13	        this->val = val;
    14	        this->rnk = rnk;
    15	    }
    16	    bool operator<(const T &b) const { return this->val < b.val; }
    17	    bool operator>(const T &b) const {
    18	        if (this->val != b.val)
    19	            return this->val > b.val;
    20	        return this->rnk < b.rnk;
    21	    }
    22	};
    23	#define param T, greater<T>, pairing_heap_tag
    24	using namespace std;
    25	using namespace __gnu_pbds;
    26	const int maxn = 1e5 + 7;
    27	__gnu_pbds::priority_queue<param> qs[maxn];
    28	int a[maxn];

    29	struct DisjointSet {
    30	    int p[maxn];
    31	    void init(int n) {
    32	        for (int i = 1; i <= n; i++) {
    33	            p[i] = i;
    34	        }
    35	    }
    36	    int find(int u) { return u == p[u] ? u : p[u] = find(p[u]); }
    37	    bool isJoined(int u, int v) { return find(u) == find(v); }

    38	    void Union(int u, int v) {
    39	        u = find(u);
    40	        v = find(v);
    41	        p[v] = u;
    42	    }
    43	} s;

    44	int main() {
    45	    int n, m;
    46	    cin >> n >> m;
    47	    s.init(n);

    48	    for (int i = 1; i <= n; i++) {
    49	        int val;
    50	        cin >> val;
    51	        qs[i].push(T(val, i));
    52	        a[i] = val;
    53	    }
    54	    for (int i = 1; i <= m; i++) {
    55	        int cmd;
    56	        cin >> cmd;
    57	        if (cmd == 1) {
    58	            int x, y;
    59	            cin >> x >> y;
    60	            if (a[x] == -1 || a[y] == -1)
    61	                continue;
    62	            if (x == y)
    63	                continue;
    64	            qs[s.find(x)].join(qs[s.find(y)]);
    65	            s.Union(x, y);
    66	        } else {
    67	            int x;
    68	            cin >> x;
    69	            if (a[x] == -1)
    70	                puts("-1");
    71	            else {
    72	                int y = s.find(x);
    73	                T tmp = qs[y].top();
    74	                cout << tmp.val << endl;
    75	                a[tmp.rnk] = -1;
    76	                qs[y].pop();
    77	            }
    78	        }
    79	    }
    80	    return 0;
    81	}
 
```
 
```cpp
//./manacher/luogu4555.cpp
     1	#include <bits/stdc++.h>
     2	#define For(i, a, b, c) for (int(i) = (a); (i) <= (b); (i) += (c))
     3	#define Bor(i, a, b, c) for (int(i) = (b); (i) >= (a); (i) -= (c))
     4	#define Min(a, b) ((a) > (b) ? (b) : (a))
     5	#define Max(a, b) ((a) > (b) ? (a) : (b))
     6	using namespace std;
     7	const int N = 200050;
     8	int p[N], ll[N], ans, rr[N], mx, id, cnt;
     9	char s[N], t[N];
    10	// 最长双回文串长度 双回文串就是两个回文串拼起来
    11	int main() {
    12	    scanf("%s", t);
    13	    int len = strlen(t);
    14	    s[++cnt] = '$', s[++cnt] = '#';
    15	    For(i, 0, len - 1, 1) s[++cnt] = t[i], s[++cnt] = '#';
    16	    s[++cnt] = '\0';
    17	    For(i, 1, cnt, 1) {
    18	        if (i < mx)
    19	            p[i] = Min(p[id * 2 - i], mx - i);
    20	        else
    21	            p[i] = 1;
    22	        while (s[i - p[i]] == s[i + p[i]])
    23	            p[i]++;
    24	        if (mx < i + p[i]) id = i, mx = i + p[i];
    25	        ll[i + p[i] - 1] = Max(ll[i + p[i] - 1], p[i] - 1);
    26	        rr[i - p[i] + 1] = Max(rr[i - p[i] + 1], p[i] - 1);
    27	    }
    28	    For(i, 2, cnt, 2) rr[i] = Max(rr[i], rr[i - 2] - 2);
    29	    Bor(i, 2, cnt, 2) ll[i] = Max(ll[i], ll[i + 2] - 2);
    30	    For(i, 2, cnt, 2) if (rr[i] && ll[i]) ans = Max(ans, ll[i] + rr[i]);
    31	    cout << ans;
    32	    return 0;
    33	}
 
```
 
```cpp
//./manacher/luogu3805.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	const int maxn = 11000002;
     4	char data[maxn << 1];
     5	int p[maxn << 1], cnt, ans;
     6	inline void qr() {
     7	    char c = getchar();
     8	    data[0] = '~', data[cnt = 1] = '|';
     9	    while (c < 'a' || c > 'z')
    10	        c = getchar();
    11	    while (c >= 'a' && c <= 'z')
    12	        data[++cnt] = c, data[++cnt] = '|', c = getchar();
    13	}
    14	// 最长回文子串 输出长度
    15	int main() {
    16	    qr();
    17	    for (int t = 1, r = 0, mid = 0; t <= cnt; ++t) {
    18	        if (t <= r) p[t] = min(p[(mid << 1) - t], r - t + 1);
    19	        while (data[t - p[t]] == data[t + p[t]])
    20	            ++p[t];
    21	        if (p[t] + t > r) r = p[t] + t - 1, mid = t;
    22	        if (p[t] > ans) ans = p[t];
    23	    }
    24	    printf("%d\n", ans - 1);
    25	    return 0;
    26	}

 
```
 
```cpp
//./AC_automaton/luogu3808.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	const int maxn = 1e6 + 6;

     4	namespace AC {
     5	const int len = 26; //字符集个数
     6	int tr[maxn][len], tot;
     7	int e[maxn], fail[maxn]; // e 代表是否是结尾
     8	void insert(string s) {
     9	    int u = 0;
    10	    for (int i = 0; i < s.length(); i++) {
    11	        int &v = tr[u][s[i] - 'a'];
    12	        if (!v) v = ++tot;
    13	        u = v;
    14	    }
    15	    e[u]++;
    16	}
    17	void build() {
    18	    queue<int> q;
    19	    for (int i = 0; i < len; i++)
    20	        if (tr[0][i]) q.push(tr[0][i]);
    21	    while (q.size()) {
    22	        int u = q.front();
    23	        q.pop();
    24	        for (int i = 0; i < len; i++) {
    25	            int &v = tr[u][i];
    26	            if (v)
    27	                fail[v] = tr[fail[u]][i], q.push(v);
    28	            else
    29	                v = tr[fail[u]][i];
    30	        }
    31	    }
    32	}
    33	int query(string t) {
    34	    int u = 0, res = 0;
    35	    for (int i = 0; i < t.length(); i++) {
    36	        u = tr[u][t[i] - 'a']; // 转移
    37	        for (int j = u; j && e[j] != -1; j = fail[j]) {
    38	            res += e[j], e[j] = -1;
    39	        }
    40	    }
    41	    return res;
    42	}
    43	} // namespace AC

    44	// 给maxn个模式串和文本串，问有多少模式串在文本串中出现

    45	int main() {
    46	    int n;
    47	    cin >> n;
    48	    string s;
    49	    for (int i = 1; i <= n; i++)
    50	        cin >> s, AC::insert(s);
    51	    cin >> s;
    52	    AC::build();
    53	    cout << AC::query(s);
    54	    return 0;
    55	}
 
```
 
```cpp
//./AC_automaton/luogu3796.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	const int N = 156, L = 1e6 + 6;
     4	namespace AC {
     5	const int SZ = N * 80;
     6	int tot, tr[SZ][26];
     7	int fail[SZ], idx[SZ], val[SZ];
     8	int cnt[N]; // 记录第 i 个字符串的出现次数
     9	void init() {
    10	    memset(fail, 0, sizeof(fail));
    11	    memset(tr, 0, sizeof(tr));
    12	    memset(val, 0, sizeof(val));
    13	    memset(cnt, 0, sizeof(cnt));
    14	    memset(idx, 0, sizeof(idx));
    15	    tot = 0;
    16	}
    17	void insert(char *s, int id) { // id 表示原始字符串的编号
    18	    int u = 0;
    19	    for (int i = 1; s[i]; i++) {
    20	        if (!tr[u][s[i] - 'a']) tr[u][s[i] - 'a'] = ++tot;
    21	        u = tr[u][s[i] - 'a'];
    22	    }
    23	    idx[u] = id;
    24	}
    25	queue<int> q;
    26	void build() {
    27	    for (int i = 0; i < 26; i++)
    28	        if (tr[0][i]) q.push(tr[0][i]);
    29	    while (q.size()) {
    30	        int u = q.front();
    31	        q.pop();
    32	        for (int i = 0; i < 26; i++) {
    33	            int v = tr[u][i];
    34	            if (v)
    35	                fail[v] = tr[fail[u]][i], q.push(v);
    36	            else
    37	                v = tr[fail[u]][i];
    38	        }
    39	    }
    40	}
    41	int query(char *t) { // 返回最大的出现次数
    42	    int u = 0, res = 0;
    43	    for (int i = 1; t[i]; i++) {
    44	        u = tr[u][t[i] - 'a'];
    45	        for (int j = u; j; j = fail[j])
    46	            val[j]++;
    47	    }
    48	    for (int i = 0; i <= tot; i++)
    49	        if (idx[i]) res = max(res, val[i]), cnt[idx[i]] = val[i];
    50	    return res;
    51	}
    52	} // namespace AC
    53	// 有 N 个由小写字母组成的模式串以及一个文本串T。
    54	// 每个模式串可能会在文本串中出现多次。你需要找出哪些模式串在文本串
    55	// T 中出现的次数最多。
    56	int n;
    57	char s[N][100], t[L];
    58	int main() {
    59	    while (~scanf("%d", &n)) {
    60	        if (n == 0) break;
    61	        AC::init();
    62	        for (int i = 1; i <= n; i++)
    63	            scanf("%s", s[i] + 1), AC::insert(s[i], i);
    64	        AC::build();
    65	        scanf("%s", t + 1);
    66	        int x = AC::query(t);
    67	        printf("%d\n", x);
    68	        for (int i = 1; i <= n; i++)
    69	            if (AC::cnt[i] == x) printf("%s\n", s[i] + 1);
    70	    }
    71	    return 0;
    72	}
 
```
 
```cpp
//./AC_automaton/POI2000.cpp
     1	#include <bits/stdc++.h>
     2	// #include <dbg_func>
     3	using namespace std;
     4	const int maxn = 1e6 + 6;
     5	namespace AC {
     6	const int len = 2; //字符集个数
     7	int tr[maxn][len], tot;
     8	int e[maxn], fail[maxn];
     9	int vis[maxn], ins[maxn];
    10	int deg[maxn], nxt[maxn], dp[maxn];
    11	char oi[maxn];
    12	void insert(string s) {
    13	    int u = 0;
    14	    for (int i = 0; i < s.length(); i++) {
    15	        int &v = tr[u][s[i] - '0'];
    16	        if (!v) v = ++tot;
    17	        oi[tot] = s[i];
    18	        u = v;
    19	    }
    20	    e[u]++;
    21	}
    22	void build() {
    23	    queue<int> q;
    24	    for (int i = 0; i < len; i++)
    25	        if (tr[0][i]) q.push(tr[0][i]);
    26	    while (q.size()) {
    27	        int u = q.front();
    28	        q.pop();
    29	        for (int i = 0; i < len; i++) {
    30	            int &v = tr[u][i];
    31	            if (v)
    32	                fail[v] = tr[fail[u]][i], e[v] |= e[fail[v]], e[v] |= e[u],
    33	                q.push(v);
    34	            else
    35	                v = tr[fail[u]][i];
    36	            // dbg(e);
    37	        }
    38	    }
    39	    if (!tr[0][0] || !tr[0][1]) {
    40	        puts("TAK");
    41	        exit(0);
    42	    }
    43	}
    44	void dfs(int now) {
    45	    if (ins[now]) {
    46	        puts("TAK");
    47	        exit(0);
    48	    }
    49	    if (vis[now] || e[now]) return;
    50	    vis[now] = ins[now] = 1;
    51	    for (int i = 0; i < len; i++)
    52	        if (tr[now][i]) {
    53	            deg[tr[now][i]]++;
    54	            dfs(tr[now][i]);
    55	        }
    56	    ins[now] = 0;
    57	}
    58	} // namespace AC

    59	// 给maxn个模式串和文本串，问有多少模式串在文本串中出现
    60	string inp[maxn];
    61	int main() {
    62	    int n;
    63	    cin >> n;
    64	    string s;
    65	    for (int i = 1; i <= n; i++) {
    66	        cin >> s;
    67	        AC::insert(s);
    68	    }
    69	    AC::build();
    70	    AC::dfs(0);
    71	    puts("NIE");
    72	    return 0;
    73	}
 
```
 
```cpp
//./AC_automaton/2016icpcHongKongJ.cpp
     1	#include <bits/stdc++.h>
     2	// #include <dbg_func>
     3	using namespace std;
     4	const int maxn = 1e6 + 6;
     5	vector<pair<int, int>> G[maxn], R[maxn];
     6	namespace AC {
     7	const int len = 2; //字符集个数
     8	int tr[maxn][len], tot;
     9	int e[maxn], fail[maxn];
    10	int vis[maxn], ins[maxn];
    11	int deg[maxn], ord[maxn], dp[maxn];
    12	int mxd[maxn];
    13	char oi[maxn];
    14	vector<int> pi;
    15	void insert(string s) {
    16	    int u = 0;
    17	    for (int i = 0; i < s.length(); i++) {
    18	        int &v = tr[u][s[i] - '0'];
    19	        if (!v) v = ++tot;
    20	        oi[tot] = s[i];
    21	        u = v;
    22	    }
    23	    e[u]++;
    24	}
    25	void build() {
    26	    queue<int> q;
    27	    for (int i = 0; i < len; i++)
    28	        if (tr[0][i]) q.push(tr[0][i]);
    29	    while (q.size()) {
    30	        int u = q.front();
    31	        q.pop();
    32	        for (int i = 0; i < len; i++) {
    33	            int &v = tr[u][i];
    34	            if (v)
    35	                fail[v] = tr[fail[u]][i], e[v] |= e[fail[v]], q.push(v);
    36	            else
    37	                v = tr[fail[u]][i];
    38	        }
    39	    }
    40	}
    41	void dfs(int now) {
    42	    vis[now] = 1;
    43	    pi.push_back(now);
    44	    for (int i = 0; i < len; i++) {
    45	        int v = tr[now][i];
    46	        if (!vis[v] && !e[v]) dfs(v);
    47	    }
    48	}

    49	int topo() {
    50	    tot = 0;
    51	    for (auto u : pi) {
    52	        for (int i = 0; i < len; i++) {
    53	            int v = tr[u][i];
    54	            if (vis[v]) {
    55	                G[u].push_back({i, v}), deg[v]++;
    56	            }
    57	        }
    58	    }
    59	    queue<int> q;
    60	    for (auto it : pi)
    61	        if (!deg[it]) q.push(it);
    62	    while (!q.empty()) {
    63	        int u = q.front();
    64	        q.pop();
    65	        ord[++tot] = u;
    66	        for (auto it : G[u]) {
    67	            if (--deg[it.second] == 0) q.push(it.second);
    68	        }
    69	    }
    70	    return tot == pi.size();
    71	}
    72	void dfs2(int u) {
    73	    mxd[u] = dp[u];
    74	    for (auto it : R[u]) {
    75	        dfs2(it.second);
    76	        mxd[u] = max(mxd[u], mxd[it.second]);
    77	    }
    78	}
    79	void print(int u) {
    80	    for (auto it : R[u]) {
    81	        if (mxd[it.second] == mxd[u]) {
    82	            cout << it.first;
    83	            print(it.second);
    84	            break;
    85	        }
    86	    }
    87	}
    88	void solve() {
    89	    for (int i = 1; i <= tot; i++) {
    90	        int u = ord[i];
    91	        for (auto it : G[u]) {
    92	            dp[it.second] = max(dp[it.second], dp[u] + 1);
    93	        }
    94	    }
    95	    int ret = 0;
    96	    for (int i = 1; i <= tot; i++)
    97	        ret = max(ret, dp[i]);
    98	    for (int i = 1; i <= tot; i++) {
    99	        int u = ord[i];
   100	        for (auto it : G[u]) {
   101	            if (dp[it.second] == dp[u] + 1) {
   102	                R[u].push_back(it);
   103	            }
   104	        }
   105	        sort(R[u].begin(), R[u].end());
   106	    }
   107	    dfs2(0);
   108	    print(0);
   109	}
   110	} // namespace AC

   111	// 2016 ICPC 香港 K 题 trie 图上找最长链
   112	string inp[maxn];
   113	int main() {
   114	    int n;
   115	    cin >> n;
   116	    string s;
   117	    for (int i = 1; i <= n; i++) {
   118	        cin >> s;
   119	        AC::insert(s);
   120	    }
   121	    AC::build();
   122	    AC::dfs(0);
   123	    if (!AC::topo()) {
   124	        cout << "-1\n";
   125	        return 0;
   126	    }
   127	    AC::solve();
   128	    return 0;
   129	}
 
```
 
```cpp
//./AC_automaton/NOI2011.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	const int maxn = 1e5 + 6;

     4	namespace AC {
     5	const int len = 26; //字符集个数
     6	int tr[maxn][len], tot;
     7	int fail[maxn], last[maxn], fa[maxn], siz[maxn], dfsq[maxn];
     8	vector<int> failtree[maxn];
     9	int enode[maxn]; // 第i个串对应的结束节点编号
    10	vector<int> e[maxn];
    11	vector<int> nxt[maxn];
    12	vector<pair<int, int>> que[maxn]; // first id second x
    13	int ans[maxn];
    14	int n, cnt, m; // dfs序
    15	int tre[maxn]; //树状数组

    16	void insert(string s) {
    17	    int u = 0;
    18	    for (int i = 0; i < s.length(); i++) {
    19	        if (s[i] == 'B') {
    20	            u = fa[u];
    21	        } else if (s[i] == 'P') {
    22	            e[u].push_back(++n);
    23	            enode[n] = u;
    24	        } else {
    25	            int &v = tr[u][s[i] - 'a'];
    26	            if (!v) {
    27	                v = ++tot;
    28	                nxt[u].push_back(tot);
    29	                fa[tot] = u;
    30	            }
    31	            u = v;
    32	        }
    33	    }
    34	}

    35	void build() {
    36	    queue<int> q;
    37	    for (int i = 0; i < len; i++)
    38	        if (tr[0][i]) q.push(tr[0][i]);
    39	    while (q.size()) {
    40	        int u = q.front();
    41	        q.pop();
    42	        for (int i = 0; i < len; i++) {
    43	            int &v = tr[u][i];
    44	            if (v)
    45	                fail[v] = tr[fail[u]][i], q.push(v);
    46	            else
    47	                v = tr[fail[u]][i];
    48	        }
    49	    }
    50	    for (int i = 1; i <= tot; i++)
    51	        failtree[fail[i]].push_back(i);
    52	}

    53	void dfs(int u) {
    54	    siz[u] = 1;
    55	    dfsq[u] = ++cnt;
    56	    for (auto i : failtree[u]) {
    57	        dfs(i);
    58	        siz[u] += siz[i];
    59	    }
    60	    // dbg(dfsq, siz);
    61	}

    62	void readquestion() {
    63	    cin >> m;
    64	    for (int i = 1; i <= m; i++) {
    65	        int x, y;
    66	        cin >> x >> y;
    67	        que[y].push_back({i, x});
    68	    }
    69	}

    70	void treeinsert(int x, int inc) {
    71	    for (x; x <= cnt; x += (x & -x))
    72	        tre[x] += inc;
    73	}

    74	int sum(int x) {
    75	    int res = 0;
    76	    for (x; x; x -= (x & -x))
    77	        res += tre[x];
    78	    return res;
    79	}

    80	void work(int u) {
    81	    treeinsert(dfsq[u], 1);
    82	    for (auto it : e[u]) {
    83	        for (auto jt : que[it]) {
    84	            int v = enode[jt.second];
    85	            ans[jt.first] += sum(dfsq[v] + siz[v] - 1) - sum(dfsq[v] - 1);
    86	        }
    87	    }
    88	    for (auto it : nxt[u])
    89	        work(it);
    90	    treeinsert(dfsq[u], -1);
    91	}

    92	void print() {
    93	    for (int i = 1; i <= m; i++)
    94	        cout << ans[i] << '\n';
    95	}
    96	} // namespace AC

    97	// 阿狸的打字机

    98	int main() {
    99	    string s;
   100	    cin >> s;
   101	    AC::insert(s);
   102	    AC::build();
   103	    AC::dfs(0);
   104	    AC::readquestion();
   105	    AC::work(0);
   106	    AC::print();
   107	    return 0;
   108	}
 
```
 
```cpp
//./AC_automaton/luogu5357.cpp
     1	#include <bits/stdc++.h>
     2	#define maxn 2000001
     3	using namespace std;
     4	char s[maxn], T[maxn];
     5	int n, cnt, vis[200051], ans, in[maxn], Map[maxn];
     6	struct kkk {
     7	    int son[26], fail, flag, ans;
     8	    void clear() {
     9	        memset(son, 0, sizeof(son)), fail = flag = ans = 0;
    10	    }
    11	} trie[maxn];
    12	queue<int> q;
    13	void insert(char *s, int num) {
    14	    int u = 1, len = strlen(s);
    15	    for (int i = 0; i < len; i++) {
    16	        int v = s[i] - 'a';
    17	        if (!trie[u].son[v]) trie[u].son[v] = ++cnt;
    18	        u = trie[u].son[v];
    19	    }
    20	    if (!trie[u].flag) trie[u].flag = num;
    21	    Map[num] = trie[u].flag;
    22	}
    23	void getFail() {
    24	    for (int i = 0; i < 26; i++)
    25	        trie[0].son[i] = 1;
    26	    q.push(1);
    27	    while (!q.empty()) {
    28	        int u = q.front();
    29	        q.pop();
    30	        int Fail = trie[u].fail;
    31	        for (int i = 0; i < 26; i++) {
    32	            int &v = trie[u].son[i];
    33	            if (!v) {
    34	                trie[u].son[i] = trie[Fail].son[i];
    35	                continue;
    36	            }
    37	            trie[v].fail = trie[Fail].son[i];
    38	            in[trie[v].fail]++;
    39	            q.push(v);
    40	        }
    41	    }
    42	}
    43	void topu() {
    44	    for (int i = 1; i <= cnt; i++)
    45	        if (in[i] == 0) q.push(i);
    46	    while (!q.empty()) {
    47	        int u = q.front();
    48	        q.pop();
    49	        vis[trie[u].flag] = trie[u].ans;
    50	        int &v = trie[u].fail;
    51	        in[v]--;
    52	        trie[v].ans += trie[u].ans;
    53	        if (in[v] == 0) q.push(v);
    54	    }
    55	}
    56	void query(char *s) {
    57	    int u = 1, len = strlen(s);
    58	    for (int i = 0; i < len; i++)
    59	        u = trie[u].son[s[i] - 'a'], trie[u].ans++;
    60	}

    61	// 给你一个文本串 S 和 n 个模式串 T_{1..n}​，
    62	// 请你分别求出每个模式串 T_i​ 在 S 中出现的次数。

    63	int main() {
    64	    scanf("%d", &n);
    65	    cnt = 1;
    66	    for (int i = 1; i <= n; i++) {
    67	        scanf("%s", s);
    68	        insert(s, i);
    69	    }
    70	    getFail();
    71	    scanf("%s", T);
    72	    query(T);
    73	    topu();
    74	    for (int i = 1; i <= n; i++)
    75	        printf("%d\n", vis[Map[i]]);
    76	}
 
```
 
```cpp
//./K-D_tree/luogu1429.cpp
     1	#include <algorithm>
     2	#include <cmath>
     3	#include <cstdio>
     4	#include <cstdlib>
     5	#include <cstring>
     6	using namespace std;
     7	const int maxn = 200010;
     8	int n, d[maxn], lc[maxn], rc[maxn];
     9	double ans = 2e18;
    10	struct node {
    11	    double x, y;
    12	} s[maxn];
    13	double L[maxn], R[maxn], D[maxn], U[maxn];
    14	double dist(int a, int b) {
    15	    return (s[a].x - s[b].x) * (s[a].x - s[b].x) +
    16	           (s[a].y - s[b].y) * (s[a].y - s[b].y);
    17	}
    18	bool cmp1(node a, node b) {
    19	    return a.x < b.x;
    20	}
    21	bool cmp2(node a, node b) {
    22	    return a.y < b.y;
    23	}
    24	void maintain(int x) {
    25	    L[x] = R[x] = s[x].x;
    26	    D[x] = U[x] = s[x].y;
    27	    if (lc[x])
    28	        L[x] = min(L[x], L[lc[x]]), R[x] = max(R[x], R[lc[x]]),
    29	        D[x] = min(D[x], D[lc[x]]), U[x] = max(U[x], U[lc[x]]);
    30	    if (rc[x])
    31	        L[x] = min(L[x], L[rc[x]]), R[x] = max(R[x], R[rc[x]]),
    32	        D[x] = min(D[x], D[rc[x]]), U[x] = max(U[x], U[rc[x]]);
    33	}
    34	int build(int l, int r) {
    35	    if (l >= r) return 0;
    36	    int mid = (l + r) >> 1;
    37	    double avx = 0, avy = 0, vax = 0, vay = 0; // average variance
    38	    for (int i = l; i <= r; i++)
    39	        avx += s[i].x, avy += s[i].y;
    40	    avx /= (double) (r - l + 1);
    41	    avy /= (double) (r - l + 1);
    42	    for (int i = l; i <= r; i++)
    43	        vax += (s[i].x - avx) * (s[i].x - avx),
    44	            vay += (s[i].y - avy) * (s[i].y - avy);
    45	    if (vax >= vay)
    46	        d[mid] = 1, nth_element(s + l, s + mid, s + r + 1, cmp1);
    47	    else
    48	        d[mid] = 2, nth_element(s + l, s + mid, s + r + 1, cmp2);
    49	    lc[mid] = build(l, mid - 1), rc[mid] = build(mid + 1, r);
    50	    maintain(mid);
    51	    return mid;
    52	}
    53	double f(int a, int b) {
    54	    double ret = 0;
    55	    if (L[b] > s[a].x) ret += (L[b] - s[a].x) * (L[b] - s[a].x);
    56	    if (R[b] < s[a].x) ret += (s[a].x - R[b]) * (s[a].x - R[b]);
    57	    if (D[b] > s[a].y) ret += (D[b] - s[a].y) * (D[b] - s[a].y);
    58	    if (U[b] < s[a].y) ret += (s[a].y - U[b]) * (s[a].y - U[b]);
    59	    return ret;
    60	}
    61	void query(int l, int r, int x) {
    62	    if (l > r) return;
    63	    int mid = (l + r) >> 1;
    64	    if (mid != x) ans = min(ans, dist(x, mid));
    65	    if (l == r) return;
    66	    double distl = f(x, lc[mid]), distr = f(x, rc[mid]);
    67	    if (distl < ans && distr < ans) {
    68	        if (distl < distr) {
    69	            query(l, mid - 1, x);
    70	            if (distr < ans) query(mid + 1, r, x);
    71	        } else {
    72	            query(mid + 1, r, x);
    73	            if (distl < ans) query(l, mid - 1, x);
    74	        }
    75	    } else {
    76	        if (distl < ans) query(l, mid - 1, x);
    77	        if (distr < ans) query(mid + 1, r, x);
    78	    }
    79	}
    80	// 平面最近点对 给定平面上的 $n$ 个点 $(x_i,y_i)$
    81	// 找出平面上最近两个点对之间的欧几里得距离
    82	int main() {
    83	    scanf("%d", &n);
    84	    for (int i = 1; i <= n; i++)
    85	        scanf("%lf%lf", &s[i].x, &s[i].y);
    86	    build(1, n);
    87	    for (int i = 1; i <= n; i++)
    88	        query(1, n, i);
    89	    printf("%.4lf\n", sqrt(ans));
    90	    return 0;
    91	}
 
```
 
```cpp
//./K-D_tree/CQOI2016.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	#define int long long
     4	const int maxn = 100010;
     5	int n, k;
     6	priority_queue<int, vector<int>, greater<int>> q;
     7	struct node {
     8	    int x, y;
     9	} s[maxn];
    10	bool cmp1(node a, node b) {
    11	    return a.x < b.x;
    12	}
    13	bool cmp2(node a, node b) {
    14	    return a.y < b.y;
    15	}
    16	int lc[maxn], rc[maxn], L[maxn], R[maxn], D[maxn], U[maxn];
    17	void maintain(int x) {
    18	    L[x] = R[x] = s[x].x;
    19	    D[x] = U[x] = s[x].y;
    20	    if (lc[x])
    21	        L[x] = min(L[x], L[lc[x]]), R[x] = max(R[x], R[lc[x]]),
    22	        D[x] = min(D[x], D[lc[x]]), U[x] = max(U[x], U[lc[x]]);
    23	    if (rc[x])
    24	        L[x] = min(L[x], L[rc[x]]), R[x] = max(R[x], R[rc[x]]),
    25	        D[x] = min(D[x], D[rc[x]]), U[x] = max(U[x], U[rc[x]]);
    26	}
    27	int build(int l, int r) {
    28	    if (l > r) return 0;
    29	    int mid = (l + r) >> 1;
    30	    double av1 = 0, av2 = 0, va1 = 0, va2 = 0; // average variance
    31	    for (int i = l; i <= r; i++)
    32	        av1 += s[i].x, av2 += s[i].y;
    33	    av1 /= (r - l + 1);
    34	    av2 /= (r - l + 1);
    35	    for (int i = l; i <= r; i++)
    36	        va1 += (av1 - s[i].x) * (av1 - s[i].x),
    37	            va2 += (av2 - s[i].y) * (av2 - s[i].y);
    38	    if (va1 > va2)
    39	        nth_element(s + l, s + mid, s + r + 1, cmp1);
    40	    else
    41	        nth_element(s + l, s + mid, s + r + 1, cmp2);
    42	    lc[mid] = build(l, mid - 1);
    43	    rc[mid] = build(mid + 1, r);
    44	    maintain(mid);
    45	    return mid;
    46	}
    47	int sq(int x) {
    48	    return x * x;
    49	}
    50	int dist(int a, int b) {
    51	    return max(sq(s[a].x - L[b]), sq(s[a].x - R[b])) +
    52	           max(sq(s[a].y - D[b]), sq(s[a].y - U[b]));
    53	}
    54	void query(int l, int r, int x) {
    55	    if (l > r) return;
    56	    int mid = (l + r) >> 1, t = sq(s[mid].x - s[x].x) + sq(s[mid].y - s[x].y);
    57	    if (t > q.top()) q.pop(), q.push(t);
    58	    int distl = dist(x, lc[mid]), distr = dist(x, rc[mid]);
    59	    if (distl > q.top() && distr > q.top()) {
    60	        if (distl > distr) {
    61	            query(l, mid - 1, x);
    62	            if (distr > q.top()) query(mid + 1, r, x);
    63	        } else {
    64	            query(mid + 1, r, x);
    65	            if (distl > q.top()) query(l, mid - 1, x);
    66	        }
    67	    } else {
    68	        if (distl > q.top()) query(l, mid - 1, x);
    69	        if (distr > q.top()) query(mid + 1, r, x);
    70	    }
    71	}
    72	// 给定平面上的 $n$ 个点 $(x_i,y_i)$ ，求欧几里得距离下的第 $k$
    73	// 远无序点对之间的距离。

    74	signed main(void) {
    75	    scanf("%lld%lld", &n, &k);
    76	    k *= 2;
    77	    for (int i = 1; i <= k; i++)
    78	        q.push(0);
    79	    for (int i = 1; i <= n; i++)
    80	        scanf("%lld%lld", &s[i].x, &s[i].y);
    81	    build(1, n);
    82	    for (int i = 1; i <= n; i++)
    83	        query(1, n, i);
    84	    printf("%lld\n", q.top());
    85	    return 0;
    86	}
 
```
 
```cpp
//./K-D_tree/luogu4148.cpp
     1	#include <algorithm>
     2	#include <cstdio>
     3	#include <cstring>
     4	using namespace std;
     5	const int maxn = 200010;
     6	int n, op, xl, xr, yl, yr, lstans;
     7	struct node {
     8	    int x, y, v;
     9	} s[maxn];
    10	bool cmp1(int a, int b) {
    11	    return s[a].x < s[b].x;
    12	}
    13	bool cmp2(int a, int b) {
    14	    return s[a].y < s[b].y;
    15	}
    16	double a = 0.725;
    17	int rt, cur, d[maxn], lc[maxn], rc[maxn], L[maxn], R[maxn], D[maxn], U[maxn],
    18	    siz[maxn], sum[maxn];
    19	int g[maxn], t;
    20	void print(int x) {
    21	    if (!x) return;
    22	    print(lc[x]);
    23	    g[++t] = x;
    24	    print(rc[x]);
    25	}
    26	void maintain(int x) {
    27	    siz[x] = siz[lc[x]] + siz[rc[x]] + 1;
    28	    sum[x] = sum[lc[x]] + sum[rc[x]] + s[x].v;
    29	    L[x] = R[x] = s[x].x;
    30	    D[x] = U[x] = s[x].y;
    31	    if (lc[x])
    32	        L[x] = min(L[x], L[lc[x]]), R[x] = max(R[x], R[lc[x]]),
    33	        D[x] = min(D[x], D[lc[x]]), U[x] = max(U[x], U[lc[x]]);
    34	    if (rc[x])
    35	        L[x] = min(L[x], L[rc[x]]), R[x] = max(R[x], R[rc[x]]),
    36	        D[x] = min(D[x], D[rc[x]]), U[x] = max(U[x], U[rc[x]]);
    37	}
    38	int build(int l, int r) {
    39	    if (l > r) return 0;
    40	    int mid = (l + r) >> 1;
    41	    double av1 = 0, av2 = 0, va1 = 0, va2 = 0;
    42	    for (int i = l; i <= r; i++)
    43	        av1 += s[g[i]].x, av2 += s[g[i]].y;
    44	    av1 /= (r - l + 1);
    45	    av2 /= (r - l + 1);
    46	    for (int i = l; i <= r; i++)
    47	        va1 += (av1 - s[g[i]].x) * (av1 - s[g[i]].x),
    48	            va2 += (av2 - s[g[i]].y) * (av2 - s[g[i]].y);
    49	    if (va1 > va2)
    50	        nth_element(g + l, g + mid, g + r + 1, cmp1), d[g[mid]] = 1;
    51	    else
    52	        nth_element(g + l, g + mid, g + r + 1, cmp2), d[g[mid]] = 2;
    53	    lc[g[mid]] = build(l, mid - 1);
    54	    rc[g[mid]] = build(mid + 1, r);
    55	    maintain(g[mid]);
    56	    return g[mid];
    57	}
    58	void rebuild(int &x) {
    59	    t = 0;
    60	    print(x);
    61	    x = build(1, t);
    62	}
    63	bool bad(int x) {
    64	    return a * siz[x] <= (double) max(siz[lc[x]], siz[rc[x]]);
    65	}
    66	void insert(int &x, int v) {
    67	    if (!x) {
    68	        x = v;
    69	        maintain(x);
    70	        return;
    71	    }
    72	    if (d[x] == 1) {
    73	        if (s[v].x <= s[x].x)
    74	            insert(lc[x], v);
    75	        else
    76	            insert(rc[x], v);
    77	    } else {
    78	        if (s[v].y <= s[x].y)
    79	            insert(lc[x], v);
    80	        else
    81	            insert(rc[x], v);
    82	    }
    83	    maintain(x);
    84	    if (bad(x)) rebuild(x);
    85	}
    86	int query(int x) {
    87	    if (!x || xr < L[x] || xl > R[x] || yr < D[x] || yl > U[x]) return 0;
    88	    if (xl <= L[x] && R[x] <= xr && yl <= D[x] && U[x] <= yr) return sum[x];
    89	    int ret = 0;
    90	    if (xl <= s[x].x && s[x].x <= xr && yl <= s[x].y && s[x].y <= yr)
    91	        ret += s[x].v;
    92	    return query(lc[x]) + query(rc[x]) + ret;
    93	}
    94	int main() {
    95	    scanf("%d", &n);
    96	    while (~scanf("%d", &op)) {
    97	        if (op == 1) {
    98	            cur++, scanf("%d%d%d", &s[cur].x, &s[cur].y, &s[cur].v);
    99	            s[cur].x ^= lstans;
   100	            s[cur].y ^= lstans;
   101	            s[cur].v ^= lstans;
   102	            insert(rt, cur);
   103	        }
   104	        if (op == 2) {
   105	            scanf("%d%d%d%d", &xl, &yl, &xr, &yr);
   106	            xl ^= lstans;
   107	            yl ^= lstans;
   108	            xr ^= lstans;
   109	            yr ^= lstans;
   110	            printf("%d\n", lstans = query(rt));
   111	        }
   112	        if (op == 3) return 0;
   113	    }
   114	}
 
```
 
```cpp
//./SA/normalSA.cpp
     1	#include <bits/stdc++.h>
     2	const int MAXN = 1e6 + 10;
     3	using namespace std;
     4	char s[MAXN];
     5	int N, M, rak[MAXN], sa[MAXN], tax[MAXN], tp[MAXN];
     6	void Qsort() {
     7	    for (int i = 0; i <= M; i++)
     8	        tax[i] = 0;
     9	    for (int i = 1; i <= N; i++)
    10	        tax[rak[i]]++;
    11	    for (int i = 1; i <= M; i++)
    12	        tax[i] += tax[i - 1];
    13	    for (int i = N; i >= 1; i--)
    14	        sa[tax[rak[tp[i]]]--] = tp[i];
    15	}
    16	void SuffixSort() {
    17	    M = 75;
    18	    for (int i = 1; i <= N; i++)
    19	        rak[i] = s[i] - '0' + 1, tp[i] = i;
    20	    Qsort();
    21	    for (int w = 1, p = 0; p < N; M = p, w <<= 1) {
    22	        p = 0;
    23	        for (int i = 1; i <= w; i++)
    24	            tp[++p] = N - w + i;
    25	        for (int i = 1; i <= N; i++)
    26	            if (sa[i] > w) tp[++p] = sa[i] - w;
    27	        Qsort();
    28	        std::swap(tp, rak);
    29	        rak[sa[1]] = p = 1;
    30	        for (int i = 2; i <= N; i++)
    31	            rak[sa[i]] = (tp[sa[i - 1]] == tp[sa[i]] &&
    32	                          tp[sa[i - 1] + w] == tp[sa[i] + w])
    33	                             ? p
    34	                             : ++p;
    35	    }
    36	    for (int i = 1; i <= N; i++)
    37	        printf("%d ", sa[i]);
    38	}
    39	int main() {
    40	    scanf("%s", s + 1);
    41	    N = strlen(s + 1);
    42	    SuffixSort();
    43	    return 0;
    44	}
 
```
 
```cpp
//./SA/luogu3809.cpp
     1	#include <assert.h>
     2	#include <stdio.h>
     3	#include <stdlib.h>
     4	#include <string.h>
     5	#define meow(args...) fprintf(stderr, args)
     6	namespace Write {
     7	const int N = 1 << 25;
     8	char buf[N], s[20], *w = buf;
     9	void flush() {
    10	    fwrite(buf, 1, w - buf, stdout);
    11	    w = buf;
    12	}
    13	inline void putchar(register char c) {
    14	    *w++ = c;
    15	}
    16	void print(register int n) {
    17	    register char *t = s;
    18	    do {
    19	        *t++ = n % 10 + 48;
    20	    } while (n /= 10);
    21	    while (t-- > s)
    22	        putchar(*t);
    23	}
    24	} // namespace Write

    25	const int N = 1e6 + 20;
    26	int ord[2 * N], sa[2 * N], rank[N], height[N], l[N], r[N], pool4[4 * N],
    27	    *mem4 = pool4;
    28	char s[N], pool1[2 * N], *mem1 = pool1;
    29	void induce(int n, int m, int *s, char *type, int *sa, int *cnt) {
    30	    int i;
    31	    memcpy(l + 1, cnt, m * sizeof(int));
    32	    memcpy(r + 1, cnt + 1, m * sizeof(int));
    33	    sa[l[s[n - 1]]++] = n - 1;
    34	    for (i = 0; i < n; ++i) {
    35	        int t = sa[i] - 1;
    36	        if (t >= 0 && type[t]) sa[l[s[t]]++] = t;
    37	    }
    38	    for (i = n; i--;) {
    39	        int t = sa[i] - 1;
    40	        if (t >= 0 && !type[t]) sa[--r[s[t]]] = t;
    41	    }
    42	}
    43	void suffix_array(int n, int m, int *s, int *sa) {
    44	    int *cnt, *lms, *s1 = s + n, *sa1 = sa + n, n1 = 0, m1 = 0, i, t;
    45	    char *type;
    46	    type = mem1;
    47	    mem1 += n + 1;
    48	    cnt = mem4;
    49	    mem4 += m + 1;
    50	    type[n] = false;
    51	    for (i = n; i--;) {
    52	        type[i] = s[i] > s[i + 1] || (s[i] == s[i + 1] && type[i + 1]);
    53	        ++cnt[s[i]];
    54	    }
    55	    for (i = 1; i <= m; ++i)
    56	        r[i] = cnt[i] += cnt[i - 1];
    57	    memset(rank, -1, n * sizeof(int));
    58	    lms = mem4;
    59	    for (i = 0; i < n; ++i) {
    60	        if (!type[i] && (i == 0 || type[i - 1])) lms[rank[i] = n1++] = i;
    61	    }
    62	    lms[n1] = n;
    63	    mem4 += n1 + 1;
    64	    memset(sa, -1, n * sizeof(int));
    65	    for (i = 0; i < n1; ++i)
    66	        sa[--r[s[lms[i]]]] = lms[i];
    67	    induce(n, m, s, type, sa, cnt);
    68	    for (i = 0, t = -1; i < n; ++i) {
    69	        int r = rank[sa[i]];
    70	        if (r != -1) {
    71	            int len = lms[r + 1] - sa[i] + 1;
    72	            m1 += t == -1 || len != lms[rank[t] + 1] - t + 1 ||
    73	                  memcmp(s + t, s + sa[i], len * sizeof(int));
    74	            s1[r] = m1;
    75	            t = sa[i];
    76	        }
    77	    }
    78	    if (n1 == m1) {
    79	        for (i = 0; i < n1; ++i)
    80	            sa1[s1[i] - 1] = i;
    81	    } else
    82	        suffix_array(n1, m1, s1, sa1);
    83	    memset(sa, -1, n * sizeof(int));
    84	    memcpy(r + 1, cnt + 1, m * sizeof(int));
    85	    for (i = n1; i--;) {
    86	        t = lms[sa1[i]];
    87	        sa[--r[s[t]]] = t;
    88	    }
    89	    induce(n, m, s, type, sa, cnt);
    90	}
    91	int main() {
    92	    int n, i, j, h;
    93	    n = fread(s, 1, N, stdin);
    94	    while (s[n - 1] > 122u)
    95	        --n;
    96	    for (i = 0; i < n; ++i) {
    97	        if (s[i] <= '9')
    98	            ord[i] = s[i] - 47;
    99	        else if (s[i] <= 'Z')
   100	            ord[i] = s[i] - 65 + 11;
   101	        else
   102	            ord[i] = s[i] - 97 + 37;
   103	    }
   104	    suffix_array(n, 26 + 26 + 10, ord, sa);
   105	    for (i = 0; i < n; ++i)
   106	        rank[sa[i]] = i;
   107	    for (i = 0, h = 0; i < n; ++i) {
   108	        if (rank[i]) {
   109	            j = sa[rank[i] - 1];
   110	            while (s[i + h] == s[j + h])
   111	                ++h;
   112	            height[rank[i]] = h;
   113	        } else
   114	            h = 0;
   115	        if (h) --h;
   116	    }
   117	    // Write::print(sa[4] + 1);
   118	    for (i = 0; i < n; ++i) {
   119	        Write::print(sa[i] + 1);
   120	        Write::putchar(" \n"[i == n - 1]);
   121	    }
   122	    Write::flush();
   123	    return 0;
   124	}
 
```
 
```cpp
//./SA/SA-IS.cpp
     1	#include <assert.h>
     2	#include <stdio.h>
     3	#include <stdlib.h>
     4	#include <string.h>
     5	#define meow(args...) fprintf(stderr, args)
     6	namespace Write {
     7	const int N = 1 << 25;
     8	char buf[N], s[20], *w = buf;
     9	void flush() {
    10	    fwrite(buf, 1, w - buf, stdout);
    11	    w = buf;
    12	}
    13	inline void putchar(register char c) {
    14	    *w++ = c;
    15	}
    16	void print(register int n) {
    17	    register char *t = s;
    18	    do {
    19	        *t++ = n % 10 + 48;
    20	    } while (n /= 10);
    21	    while (t-- > s)
    22	        putchar(*t);
    23	}
    24	} // namespace Write

    25	const int N = 1e5 + 20;
    26	int ord[2 * N], sa[2 * N], rank[N], height[N], l[N], r[N], pool4[4 * N],
    27	    *mem4 = pool4;
    28	char s[N], pool1[2 * N], *mem1 = pool1;
    29	void induce(int n, int m, int *s, char *type, int *sa, int *cnt) {
    30	    int i;
    31	    memcpy(l + 1, cnt, m * sizeof(int));
    32	    memcpy(r + 1, cnt + 1, m * sizeof(int));
    33	    sa[l[s[n - 1]]++] = n - 1;
    34	    for (i = 0; i < n; ++i) {
    35	        int t = sa[i] - 1;
    36	        if (t >= 0 && type[t]) sa[l[s[t]]++] = t;
    37	    }
    38	    for (i = n; i--;) {
    39	        int t = sa[i] - 1;
    40	        if (t >= 0 && !type[t]) sa[--r[s[t]]] = t;
    41	    }
    42	}
    43	void suffix_array(int n, int m, int *s, int *sa) {
    44	    int *cnt, *lms, *s1 = s + n, *sa1 = sa + n, n1 = 0, m1 = 0, i, t;
    45	    char *type;
    46	    type = mem1;
    47	    mem1 += n + 1;
    48	    cnt = mem4;
    49	    mem4 += m + 1;
    50	    type[n] = false;
    51	    for (i = n; i--;) {
    52	        type[i] = s[i] > s[i + 1] || (s[i] == s[i + 1] && type[i + 1]);
    53	        ++cnt[s[i]];
    54	    }
    55	    for (i = 1; i <= m; ++i)
    56	        r[i] = cnt[i] += cnt[i - 1];
    57	    memset(rank, -1, n * sizeof(int));
    58	    lms = mem4;
    59	    for (i = 0; i < n; ++i) {
    60	        if (!type[i] && (i == 0 || type[i - 1])) lms[rank[i] = n1++] = i;
    61	    }
    62	    lms[n1] = n;
    63	    mem4 += n1 + 1;
    64	    memset(sa, -1, n * sizeof(int));
    65	    for (i = 0; i < n1; ++i)
    66	        sa[--r[s[lms[i]]]] = lms[i];
    67	    induce(n, m, s, type, sa, cnt);
    68	    for (i = 0, t = -1; i < n; ++i) {
    69	        int r = rank[sa[i]];
    70	        if (r != -1) {
    71	            int len = lms[r + 1] - sa[i] + 1;
    72	            m1 += t == -1 || len != lms[rank[t] + 1] - t + 1 ||
    73	                  memcmp(s + t, s + sa[i], len * sizeof(int));
    74	            s1[r] = m1;
    75	            t = sa[i];
    76	        }
    77	    }
    78	    if (n1 == m1) {
    79	        for (i = 0; i < n1; ++i)
    80	            sa1[s1[i] - 1] = i;
    81	    } else
    82	        suffix_array(n1, m1, s1, sa1);
    83	    memset(sa, -1, n * sizeof(int));
    84	    memcpy(r + 1, cnt + 1, m * sizeof(int));
    85	    for (i = n1; i--;) {
    86	        t = lms[sa1[i]];
    87	        sa[--r[s[t]]] = t;
    88	    }
    89	    induce(n, m, s, type, sa, cnt);
    90	}

    91	// 读入一个长度为 n 的由小写英文字母组成的字符串，
    92	// 请把这个字符串的所有非空后缀按字典序从小到大排序，
    93	// 然后按顺序输出后缀的第一个字符在原串中的位置。位置编号为 1 到 n。
    94	// 除此之外为了进一步证明你确实有给后缀排序的超能力，请另外输出
    95	// n−1 个整数分别表示排序后相邻后缀的最长公共前缀的长度。

    96	int main() {
    97	    int n, i, j, h;
    98	    n = fread(s, 1, N, stdin);
    99	    while (s[n - 1] - 97u > 25)
   100	        --n;
   101	    for (i = 0; i < n; ++i)
   102	        ord[i] = s[i] - 96;
   103	    suffix_array(n, 26, ord, sa);
   104	    for (i = 0; i < n; ++i)
   105	        rank[sa[i]] = i;
   106	    for (i = 0, h = 0; i < n; ++i) {
   107	        if (rank[i]) {
   108	            j = sa[rank[i] - 1];
   109	            while (s[i + h] == s[j + h])
   110	                ++h;
   111	            height[rank[i]] = h;
   112	        } else
   113	            h = 0;
   114	        if (h) --h;
   115	    }
   116	    for (i = 0; i < n; ++i) {
   117	        Write::print(sa[i] + 1);
   118	        Write::putchar(" \n"[i == n - 1]);
   119	    }
   120	    for (i = 1; i < n; ++i) {
   121	        Write::print(height[i]);
   122	        Write::putchar(" \n"[i == n - 1]);
   123	    }
   124	    Write::flush();
   125	    return 0;
   126	}
 
```
 
```cpp
//./sequential_automaton/HEOI2015.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	const int N = 2005;
     4	char s[N], t[N];
     5	int na[N][26], nb[N][26], nxt[26];
     6	int n, m, a[N], b[N], tot = 1, p = 1, f[N][N << 1];
     7	struct SAM {
     8	    int par, ch[26], len;
     9	} sam[N << 1];
    10	void insert(int x) {
    11	    int np = ++tot;
    12	    sam[np].len = sam[p].len + 1;
    13	    while (p && !sam[p].ch[x]) {
    14	        sam[p].ch[x] = np;
    15	        p = sam[p].par;
    16	    }
    17	    if (p == 0)
    18	        sam[np].par = 1;
    19	    else {
    20	        int q = sam[p].ch[x];
    21	        if (sam[q].len == sam[p].len + 1)
    22	            sam[np].par = q;
    23	        else {
    24	            int nq = ++tot;
    25	            sam[nq].len = sam[p].len + 1;
    26	            memcpy(sam[nq].ch, sam[q].ch, sizeof(sam[q].ch));
    27	            sam[nq].par = sam[q].par;
    28	            sam[q].par = sam[np].par = nq;

    29	            while (p && sam[p].ch[x] == q) {
    30	                sam[p].ch[x] = nq;
    31	                p = sam[p].par;
    32	            }
    33	        }
    34	    }
    35	    p = np;
    36	}
    37	int main() {
    38	    scanf("%s%s", s + 1, t + 1);
    39	    n = strlen(s + 1);
    40	    m = strlen(t + 1);
    41	    for (int i = 1; i <= n; ++i)
    42	        a[i] = s[i] - 'a';
    43	    for (int i = 1; i <= m; ++i)
    44	        b[i] = t[i] - 'a';
    45	    for (int i = 1; i <= m; ++i)
    46	        insert(b[i]);
    47	    for (int i = 0; i < 26; ++i)
    48	        nxt[i] = n + 1;
    49	    for (int i = n; i >= 0; --i) {
    50	        memcpy(na[i], nxt, sizeof(nxt));
    51	        nxt[a[i]] = i;
    52	    }
    53	    for (int i = 0; i < 26; ++i)
    54	        nxt[i] = m + 1;
    55	    for (int i = m; i >= 0; --i) {
    56	        memcpy(nb[i], nxt, sizeof(nxt));
    57	        nxt[b[i]] = i;
    58	    }
    59	    int ans = N;
    60	    for (int l = 1; l <= n; ++l) {
    61	        for (int r = l, u = 1; r <= n; ++r) {
    62	            u = sam[u].ch[a[r]];

    63	            if (!u) {
    64	                ans = min(ans, r - l + 1);
    65	                break;
    66	            }
    67	        }
    68	    }
    69	    printf("%d\n", ans == N ? -1 : ans);
    70	    ans = N;
    71	    for (int l = 1; l <= n; ++l) {
    72	        for (int r = l, u = 0; r <= n; ++r) {
    73	            u = nb[u][a[r]];

    74	            if (u == m + 1) {
    75	                ans = min(ans, r - l + 1);
    76	                break;
    77	            }
    78	        }
    79	    }
    80	    printf("%d\n", ans == N ? -1 : ans);
    81	    for (int i = n; i >= 0; --i) {
    82	        for (int j = 1; j <= tot; ++j) {
    83	            f[i][j] = N;

    84	            for (int c = 0; c < 26; ++c) {
    85	                int u = na[i][c];
    86	                int v = sam[j].ch[c];

    87	                if (u <= n) f[i][j] = min(f[i][j], f[u][v] + 1);
    88	            }
    89	        }
    90	    }
    91	    printf("%d\n", f[0][1] == N ? -1 : f[0][1]);
    92	    memset(f, 0, sizeof(f));
    93	    for (int i = n; i >= 0; --i) {
    94	        for (int j = 0; j <= m; ++j) {
    95	            f[i][j] = N;

    96	            for (int c = 0; c < 26; ++c) {
    97	                int u = na[i][c];
    98	                int v = nb[j][c];

    99	                if (u <= n) f[i][j] = min(f[i][j], f[u][v] + 1);
   100	            }
   101	        }
   102	    }
   103	    printf("%d\n", f[0][0] == N ? -1 : f[0][0]);
   104	    return 0;
   105	}
 
```
 
```cpp
//./minimum_ball_coverage/poj2069.cpp
     1	#include <bits/stdc++.h>
     2	#define T 100
     3	#define eps 1e-8
     4	#define delta 0.98
     5	#define INF 0x7fffffff
     6	using namespace std;
     7	typedef long long ll;
     8	const int maxn = 1e2 + 5;

     9	struct point {
    10	    double x, y, z;
    11	} p[maxn];

    12	double dis(point A, point B) {
    13	    return sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
    14	                (A.z - B.z) * (A.z - B.z));
    15	}

    16	double GetSum(int n, point t) {
    17	    double ans = 0;
    18	    for (int i = 0; i < n; i++)
    19	        ans += dis(t, p[i]);
    20	    return ans;
    21	}

    22	double Search(int n) {
    23	    point s = p[0];
    24	    double t = T;
    25	    double ans = INF;
    26	    while (t > eps) {
    27	        double maxx = 0;
    28	        int k = 0;
    29	        for (int i = 0; i < n; i++) {
    30	            double temp = dis(p[i], s);
    31	            if (temp > maxx) {
    32	                maxx = temp;
    33	                k = i;
    34	            }
    35	        }
    36	        ans = min(ans, maxx);
    37	        s.x += (p[k].x - s.x) / maxx * t;
    38	        s.y += (p[k].y - s.y) / maxx * t;
    39	        s.z += (p[k].z - s.z) / maxx * t;
    40	        t *= delta;
    41	    }
    42	    return ans;
    43	}

    44	// 朴素的板子题 给你一堆点，求最小球覆盖

    45	int main() {
    46	    int n;
    47	    while (scanf("%d", &n) && n) {
    48	        for (int i = 0; i < n; i++)
    49	            scanf("%lf%lf%lf", &p[i].x, &p[i].y, &p[i].z);
    50	        double ans = Search(n);
    51	        printf("%.5f\n", ans);
    52	    }
    53	    return 0;
    54	}
 
```
 
```cpp
//./BM/bm.cpp
     1	// BM-线性递推
     2	#include <algorithm>
     3	#include <bits/stdc++.h>
     4	#include <cassert>
     5	#include <cmath>
     6	#include <cstdio>
     7	#include <cstring>
     8	#include <map>
     9	#include <set>
    10	#include <string>
    11	#include <vector>
    12	using namespace std;
    13	#define rep(i, a, n) for (int i = a; i < n; i++)
    14	#define per(i, a, n) for (int i = n - 1; i >= a; i--)
    15	#define pb push_back
    16	#define mp make_pair
    17	#define all(x) (x).begin(), (x).end()
    18	#define fi first
    19	#define se second
    20	#define SZ(x) ((int) (x).size())
    21	typedef vector<int> VI;
    22	typedef long long ll;
    23	typedef pair<int, int> PII;
    24	const ll mod = 1000000007;
    25	ll powmod(ll a, ll b) {
    26	    ll res = 1;
    27	    a %= mod;
    28	    assert(b >= 0);
    29	    for (; b; b >>= 1) {
    30	        if (b & 1) res = res * a % mod;
    31	        a = a * a % mod;
    32	    }
    33	    return res;
    34	}
    35	// head

    36	long long _, n;
    37	namespace linear_seq {
    38	const int N = 10010;
    39	ll res[N], base[N], _c[N], _md[N];

    40	vector<int> Md;
    41	void mul(ll *a, ll *b, int k) {
    42	    rep(i, 0, k + k) _c[i] = 0;
    43	    rep(i, 0, k) if (a[i]) rep(j, 0, k) _c[i + j] =
    44	        (_c[i + j] + a[i] * b[j]) % mod;
    45	    for (int i = k + k - 1; i >= k; i--)
    46	        if (_c[i])
    47	            rep(j, 0, SZ(Md)) _c[i - k + Md[j]] =
    48	                (_c[i - k + Md[j]] - _c[i] * _md[Md[j]]) % mod;
    49	    rep(i, 0, k) a[i] = _c[i];
    50	}
    51	int solve(ll n, VI a, VI b) { // a 系数 b 初值 b[n+1]=a[0]*b[n]+...
    52	                              //        printf("%d\n",SZ(b));
    53	    ll ans = 0, pnt = 0;
    54	    int k = SZ(a);
    55	    assert(SZ(a) == SZ(b));
    56	    rep(i, 0, k) _md[k - 1 - i] = -a[i];
    57	    _md[k] = 1;
    58	    Md.clear();
    59	    rep(i, 0, k) if (_md[i] != 0) Md.push_back(i);
    60	    rep(i, 0, k) res[i] = base[i] = 0;
    61	    res[0] = 1;
    62	    while ((1ll << pnt) <= n)
    63	        pnt++;
    64	    for (int p = pnt; p >= 0; p--) {
    65	        mul(res, res, k);
    66	        if ((n >> p) & 1) {
    67	            for (int i = k - 1; i >= 0; i--)
    68	                res[i + 1] = res[i];
    69	            res[0] = 0;
    70	            rep(j, 0, SZ(Md)) res[Md[j]] =
    71	                (res[Md[j]] - res[k] * _md[Md[j]]) % mod;
    72	        }
    73	    }
    74	    rep(i, 0, k) ans = (ans + res[i] * b[i]) % mod;
    75	    if (ans < 0) ans += mod;
    76	    return ans;
    77	}
    78	VI BM(VI s) {
    79	    VI C(1, 1), B(1, 1);
    80	    int L = 0, m = 1, b = 1;
    81	    rep(n, 0, SZ(s)) {
    82	        ll d = 0;
    83	        rep(i, 0, L + 1) d = (d + (ll) C[i] * s[n - i]) % mod;
    84	        if (d == 0)
    85	            ++m;
    86	        else if (2 * L <= n) {
    87	            VI T = C;
    88	            ll c = mod - d * powmod(b, mod - 2) % mod;
    89	            while (SZ(C) < SZ(B) + m)
    90	                C.pb(0);
    91	            rep(i, 0, SZ(B)) C[i + m] = (C[i + m] + c * B[i]) % mod;
    92	            L = n + 1 - L;
    93	            B = T;
    94	            b = d;
    95	            m = 1;
    96	        } else {
    97	            ll c = mod - d * powmod(b, mod - 2) % mod;
    98	            while (SZ(C) < SZ(B) + m)
    99	                C.pb(0);
   100	            rep(i, 0, SZ(B)) C[i + m] = (C[i + m] + c * B[i]) % mod;
   101	            ++m;
   102	        }
   103	    }
   104	    return C;
   105	}
   106	int gao(VI a, ll n) {
   107	    VI c = BM(a);
   108	    c.erase(c.begin());
   109	    rep(i, 0, SZ(c)) c[i] = (mod - c[i]) % mod;
   110	    return solve(n, c, VI(a.begin(), a.begin() + SZ(c)));
   111	}
   112	}; // namespace linear_seq

   113	int main() {
   114	    while (~scanf("%lld", &n)) {
   115	        // 在这里输入前n项
   116	        printf("%d\n",
   117	               linear_seq::gao(
   118	                   vector<int>{1,       1,       5,        11,       36,
   119	                               95,      281,     781,      2245,     6336,
   120	                               18061,   51205,   145601,   413351,   1174500,
   121	                               3335651, 9475901, 26915305, 76455961, 217172736},
   122	                   n));
   123	    }
   124	}
 
```
 
```cpp
//./ex_kmp/hdu6192.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	static auto fast_io = []() {
     4	    std::ios::sync_with_stdio(false); // turn off sync
     5	    cin.tie(nullptr);                 // untie in/out streams
     6	    return 0;
     7	}();
     8	vector<int> exkmp(string s) {
     9	    int n = (int) s.length();
    10	    vector<int> z(n);
    11	    for (int i = 1, l = 0, r = 0; i < n; ++i) {
    12	        if (i <= r) z[i] = min(r - i + 1, z[i - l]);
    13	        while (i + z[i] < n && s[z[i]] == s[i + z[i]])
    14	            ++z[i];
    15	        if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
    16	    }
    17	    return z;
    18	}
    19	const int maxn = 1e5 + 7;
    20	const unsigned int mi = 131;
    21	unsigned int ex[maxn], has[maxn];
    22	string s;
    23	unsigned int _seg(int l, int r, int k, char x) {
    24	    return has[r] - (l - 1 >= 0 ? has[l - 1] : 0) * ex[r - l + 1] -
    25	           (k >= l && k <= r ? ex[r - k] * s[k] - ex[r - k] * x : 0);
    26	};
    27	// 修改字符串中一个字符，使得循环节最短（最后一个可以不取完全）
    28	int main() {
    29	    int n;
    30	    while (cin >> n >> s) {
    31	        // dbg(n, s);
    32	        vector<int> s_exkmp = exkmp(s);
    33	        ex[0] = 1, has[0] = s[0];
    34	        // dbg(s_exkmp);
    35	        for (int i = 1; i <= n; i++)
    36	            ex[i] = ex[i - 1] * mi;
    37	        for (int i = 1; i < n; i++)
    38	            has[i] = has[i - 1] * mi + s[i];
    39	        for (int i = 1; i < n; i++) {
    40	            int v = s_exkmp[i];
    41	            if (v + i == n) {
    42	                cout << i << ' ' << n << '\n';
    43	                break;
    44	            } else {
    45	                int p = 0;
    46	                // dbg(i);
    47	                p += (_seg(0, n - 1 - i, v, s[i + v]) ==
    48	                      _seg(i, n - 1, v, s[i + v]));
    49	                p += (_seg(0, n - 1 - i, i + v, s[v]) ==
    50	                      _seg(i, n - 1, i + v, s[v]));
    51	                if (p) {
    52	                    cout << i << ' ' << p << '\n';
    53	                    break;
    54	                }
    55	            }
    56	        }
    57	        if (n == 1) cout << "1 1\n";
    58	    }
    59	}
 
```
 
```cpp
//./ex_kmp/luogu5410.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	static auto fast_io = []() {
     4	    std::ios::sync_with_stdio(false); // turn off sync
     5	    cin.tie(nullptr);                 // untie in/out streams
     6	    return 0;
     7	}();
     8	vector<int> exkmp(string s) {
     9	    int n = (int) s.length();
    10	    vector<int> z(n);
    11	    for (int i = 1, l = 0, r = 0; i < n; ++i) {
    12	        if (i <= r) z[i] = min(r - i + 1, z[i - l]);
    13	        while (i + z[i] < n && s[z[i]] == s[i + z[i]])
    14	            ++z[i];
    15	        if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
    16	    }
    17	    return z;
    18	}
    19	// 给两个字符串a,b，求b的exkmp和b与a的每个后缀的LCP
    20	int main() {
    21	    string a, b;
    22	    cin >> a >> b;
    23	    string c = b + '?' + a;
    24	    vector<int> b_exkmp = exkmp(b);
    25	    vector<int> c_exkmp = exkmp(c);
    26	    int lenb = b.length(), lena = a.length();
    27	    long long ans1 = lenb + 1, ans2 = c_exkmp[lenb + 1] + 1;
    28	    for (int i = 1; i < lenb; i++) {
    29	        ans1 ^= (i + 1LL) * (b_exkmp[i] + 1LL);
    30	    }
    31	    for (int i = 2; i <= lena; i++) {
    32	        ans2 ^= i * (c_exkmp[i + lenb] + 1LL);
    33	    }
    34	    cout << ans1 << '\n' << ans2;
    35	}
 
```
 
```cpp
//./mo_algorithm/luogu1494.cpp
     1	#include <algorithm>
     2	#include <cmath>
     3	#include <cstdio>
     4	using namespace std;
     5	const int N = 50005;
     6	int n, m, maxn;
     7	int c[N];
     8	long long sum;
     9	int cnt[N];
    10	long long ans1[N], ans2[N];
    11	struct query {
    12	    int l, r, id;
    13	    bool operator<(const query &x) const {
    14	        if (l / maxn != x.l / maxn) return l < x.l;
    15	        return (l / maxn) & 1 ? r < x.r : r > x.r;
    16	    }
    17	} a[N];
    18	void add(int i) {
    19	    sum += cnt[i];
    20	    cnt[i]++;
    21	}
    22	void del(int i) {
    23	    cnt[i]--;
    24	    sum -= cnt[i];
    25	}
    26	long long gcd(long long a, long long b) {
    27	    return b ? gcd(b, a % b) : a;
    28	}
    29	int main() {
    30	    scanf("%d%d", &n, &m);
    31	    maxn = sqrt(n);
    32	    for (int i = 1; i <= n; i++)
    33	        scanf("%d", &c[i]);
    34	    for (int i = 0; i < m; i++)
    35	        scanf("%d%d", &a[i].l, &a[i].r), a[i].id = i;
    36	    sort(a, a + m);
    37	    for (int i = 0, l = 1, r = 0; i < m; i++) {
    38	        if (a[i].l == a[i].r) {
    39	            ans1[a[i].id] = 0, ans2[a[i].id] = 1;
    40	            continue;
    41	        }
    42	        while (l > a[i].l)
    43	            add(c[--l]);
    44	        while (r < a[i].r)
    45	            add(c[++r]);
    46	        while (l < a[i].l)
    47	            del(c[l++]);
    48	        while (r > a[i].r)
    49	            del(c[r--]);
    50	        ans1[a[i].id] = sum;
    51	        ans2[a[i].id] = (long long) (r - l + 1) * (r - l) / 2;
    52	    }
    53	    for (int i = 0; i < m; i++) {
    54	        if (ans1[i] != 0) {
    55	            long long g = gcd(ans1[i], ans2[i]);
    56	            ans1[i] /= g, ans2[i] /= g;
    57	        } else
    58	            ans2[i] = 1;
    59	        printf("%lld/%lld\n", ans1[i], ans2[i]);
    60	    }
    61	    return 0;
    62	}

 
```
 
```cpp
//./suffix_automaton/luogu3804.cpp
     1	#include <bits/stdc++.h>
     2	#ifndef open_dbg_func
     3	#define dbg(args...) (args)
     4	#endif
     5	using namespace std;
     6	struct state {
     7	    int len, link;
     8	    map<char, int> nxt;
     9	};
    10	const int maxn = 1e6 + 7;
    11	long long ans;
    12	state st[maxn * 2];
    13	int len[maxn * 2], siz[maxn * 2];
    14	int sz, last;
    15	vector<int> leni[maxn];
    16	void sam_init() {
    17	    for (int i = 0; i < sz; i++) {
    18	        st[i].nxt.clear();
    19	        siz[i] = 0;
    20	    }
    21	    st[0].len = 0;
    22	    st[0].link = -1;
    23	    sz = 1;
    24	    last = 0;
    25	}

    26	void sam_extend(char c) {
    27	    int cur = sz++;
    28	    st[cur].len = st[last].len + 1;
    29	    int p = last;
    30	    while (p != -1 && !st[p].nxt.count(c)) {
    31	        st[p].nxt[c] = cur;
    32	        p = st[p].link;
    33	    }
    34	    siz[cur] = 1;
    35	    if (p == -1) {
    36	        st[cur].link = 0;
    37	    } else {
    38	        int q = st[p].nxt[c];
    39	        if (st[p].len + 1 == st[q].len) {
    40	            st[cur].link = q;
    41	        } else {
    42	            int clone = sz++;
    43	            st[clone].len = st[p].len + 1;
    44	            st[clone].nxt = st[q].nxt;
    45	            st[clone].link = st[q].link;
    46	            while (p != -1 && st[p].nxt[c] == q) {
    47	                st[p].nxt[c] = clone;
    48	                p = st[p].link;
    49	            }
    50	            st[q].link = st[cur].link = clone;
    51	        }
    52	    }
    53	    last = cur;
    54	}
    55	int main() {
    56	    string s;
    57	    cin >> s;
    58	    sam_init();
    59	    for (int i = 0; i < s.length(); i++)
    60	        sam_extend(s[i]);
    61	    for (int i = 1; i < sz; i++) {
    62	        leni[st[i].len].push_back(i);
    63	    }
    64	    ans = 0;
    65	    for (int i = s.length(); i > 0; i--) {
    66	        for (auto it : leni[i]) {
    67	            siz[st[it].link] += siz[it];
    68	            if (siz[it] > 1) ans = max(ans, 1LL * siz[it] * st[it].len);
    69	        }
    70	    }
    71	    cout << ans;
    72	}
 
```
 
```cpp
//./suffix_automaton/std.cpp
     1	#include <bits/stdc++.h>
     2	#define N 2000005
     3	typedef long long ll;
     4	using namespace std;
     5	char s[N];
     6	int a[N], c[N], size[N], n;
     7	ll ans = 0;
     8	struct SuffixAutoMaton {
     9	    int last, cnt;
    10	    int ch[N << 1][26], fa[N << 1], l[N << 1];
    11	    void ins(int c) {
    12	        int p = last, np = ++cnt;
    13	        last = np;
    14	        l[np] = l[p] + 1;
    15	        for (; p && !ch[p][c]; p = fa[p])
    16	            ch[p][c] = np;
    17	        if (!p)
    18	            fa[np] = 1;
    19	        else {
    20	            int q = ch[p][c];
    21	            if (l[p] + 1 == l[q])
    22	                fa[np] = q;
    23	            else {
    24	                int nq = ++cnt;
    25	                l[nq] = l[p] + 1;
    26	                memcpy(ch[nq], ch[q], sizeof(ch[q]));
    27	                fa[nq] = fa[q];
    28	                fa[q] = fa[np] = nq;
    29	                for (; ch[p][c] == q; p = fa[p])
    30	                    ch[p][c] = nq;
    31	            }
    32	        }
    33	        size[np] = 1;
    34	    }
    35	    void build() {
    36	        scanf("%s", s + 1);
    37	        int len = strlen(s + 1);
    38	        last = cnt = 1;
    39	        for (int i = 1; i <= len; i++)
    40	            ins(s[i] - 'a');
    41	    }
    42	    void calc() {
    43	        for (int i = 1; i <= cnt; i++)
    44	            c[l[i]]++;
    45	        for (int i = 1; i <= cnt; i++)
    46	            c[i] += c[i - 1];
    47	        for (int i = 1; i <= cnt; i++)
    48	            a[c[l[i]]--] = i;
    49	        for (int i = cnt; i; i--) {
    50	            int p = a[i];
    51	            size[fa[p]] += size[p];
    52	            if (size[p] > 1) ans = max(ans, 1LL * size[p] * l[p]);
    53	        }
    54	        printf("%lld\n", ans);
    55	    }
    56	} sam;
    57	int main() {
    58	    sam.build();
    59	    sam.calc();
    60	}
 
```
 
```cpp
//./Gaussian_elimination/luogu3389.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	const int maxn = 505;
     4	const double eps = 1e-8;
     5	double A[maxn << 1][maxn],
     6	    x[maxn]; // A矩阵中每一行1~n存系数，n+1为答案，m个方程m行，x是最终的答案
     7	//注意空间要多开几个，还要考虑n，m不同的情况

     8	int Guass(int n, int m) //有n个未知数，m个方程
     9	{
    10	    int i = 1, j = 1, k, r, c;
    11	    while (i <= m && j <= n) //正在处理第i个方程，解第j个未知数
    12	    {
    13	        r = i; //找到绝对值最大的系数，防止除数为0的情况，使得其他方程组系数不会变得太大
    14	        for (k = i + 1; k <= m; k++)
    15	            if (fabs(A[k][j]) > fabs(A[r][j])) r = k;
    16	        if (fabs(A[r][j]) >=
    17	            eps) //出现为0的情况，说明此项已经被消掉了，直接用进行下一个未知数，而方程不变,不过这个时候，一般来说跳过的这个元素就没有固定解啦
    18	        {
    19	            for (c = 1; c <= n + 1; c++)
    20	                swap(A[i][c], A[r][c]); //交换
    21	            for (k = i + 1; k <= m; k++)
    22	                if (fabs(A[k][j]) >= eps) {
    23	                    double f = A[k][j] / A[i][j];
    24	                    for (c = j; c <= n + 1; c++) //当前方程j前面的系数都是0
    25	                        A[k][c] -= f * A[i][c];
    26	                }
    27	            i++; //获取下一个方程
    28	        }
    29	        j++; //去消下一个未知数
    30	    }
    31	    //必须先判无解再判断多解
    32	    for (k = i; k <= m; k++)
    33	        if (fabs(A[k][n + 1]) >= eps)
    34	            return 0; //若有一行系数为0但是不为答案，则无解
    35	    if (i <= n)
    36	        return 2; //如果被你处理出来的方程没有n个，就会出现多解。(i=n表示解决了n-1个方程)
    37	    for (int i = n; i >= 1; i--) {
    38	        for (j = i + 1; j <= n; j++)
    39	            A[i][n + 1] -= A[i][j] * x[j];
    40	        x[i] = A[i][n + 1] / A[i][i];
    41	    }
    42	    //最终统计出来的答案x[i]肯定是对应的第i个元素的解哦,换的只是方程的顺序
    43	    return 1; //拥有唯一解
    44	}
    45	// 模线性方程组 除法换成逆元

    46	int B[maxn << 1][maxn], cnt[maxn], ans[maxn][maxn];
    47	void xorGauss(int n, int m) {
    48	    int i = 1, j = 1, k, r, c;
    49	    while (i <= m && j <= n) // i是正在考虑的方程，j是待求解的系数
    50	    {
    51	        for (r = i; r <= m; r++)
    52	            if (B[r][j]) break; //找到有值的数简单些
    53	        if (r <= m)             //新判断
    54	        {
    55	            for (c = 1; c <= n + 1; c++)
    56	                swap(B[r][c], B[i][c]);
    57	            for (k = i + 1; k <= m; k++)
    58	                if (B[k][j])
    59	                    for (c = j; c <= n + 1; c++)
    60	                        B[k][c] ^= B[i][c];
    61	            i++;
    62	        }
    63	        j++;
    64	    }
    65	    i--, j--;
    66	    //上面时解方程部分，下面是搜索自由元部分
    67	    /*
    68	    不搜索自由元可以这样
    69	    for(k=i;k<=m;k++)if(B[k][n+1])return -1;//判断无解
    70	    return n-i+1;//返回自由元个数
    71	    */
    72	    for (int k = 1; k <= i; k++)
    73	        for (int j = n; j >= 1; j--)
    74	            if (B[k][j]) cnt[k] = j;
    75	    while (
    76	        i >= 1 ||
    77	        j >=
    78	            1) //即使没有方程了，也要把自由元算完，但是不会出现还有方程却没有自由元的情况
    79	    {
    80	        if (cnt[i] != j) {
    81	            ans[1][j] = 1;
    82	            j--;
    83	        } else {
    84	            for (int k = j + 1; k <= n; k++)
    85	                if (B[i][k]) B[i][n + 1] ^= ans[1][k];
    86	            ans[1][j] = B[i][n + 1];
    87	            i--, j--;
    88	        }
    89	    }
    90	}
    91	// 异或方程组

    92	int main() {
    93	    int n, m;
    94	    cin >> n;
    95	    m = n;
    96	    for (int i = 1; i <= m; i++)
    97	        for (int j = 1; j <= n + 1; j++)
    98	            cin >> A[i][j];
    99	    int ans = Guass(n, m);
   100	    if (ans != 1)
   101	        cout << "No Solution\n";
   102	    else
   103	        for (int i = 1; i <= n; i++)
   104	            cout << fixed << setprecision(2) << x[i] << '\n';
   105	}
 
```
 
```cpp
//./inverse/hdu4773.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;

     3	const double EPS = 1e-8;      // 精度系数
     4	const double PI = acos(-1.0); // π
     5	const int N = 4;

     6	struct Point {
     7	    double x, y;
     8	    Point(double x = 0, double y = 0)
     9	        : x(x)
    10	        , y(y) {
    11	    }
    12	    const bool operator<(Point A) const {
    13	        return x == A.x ? y < A.y : x < A.x;
    14	    }
    15	}; // 点的定义

    16	typedef Point Vector; // 向量的定义

    17	Vector operator+(Vector A, Vector B) {
    18	    return Vector(A.x + B.x, A.y + B.y);
    19	} // 向量加法
    20	Vector operator-(Vector A, Vector B) {
    21	    return Vector(A.x - B.x, A.y - B.y);
    22	} // 向量减法
    23	Vector operator*(Vector A, double p) {
    24	    return Vector(A.x * p, A.y * p);
    25	} // 向量数乘
    26	Vector operator/(Vector A, double p) {
    27	    return Vector(A.x / p, A.y / p);
    28	} // 向量数除

    29	int dcmp(double x) {
    30	    if (fabs(x) < EPS)
    31	        return 0;
    32	    else
    33	        return x < 0 ? -1 : 1;
    34	} // 与0的关系

    35	double Dot(Vector A, Vector B) {
    36	    return A.x * B.x + A.y * B.y;
    37	} // 向量点乘
    38	double Length(Vector A) {
    39	    return sqrt(Dot(A, A));
    40	} // 向量长度
    41	double Cross(Vector A, Vector B) {
    42	    return A.x * B.y - A.y * B.x;
    43	} // 向量叉乘

    44	Point GetLineProjection(Point P, Point A, Point B) {
    45	    Vector v = B - A;
    46	    return A + v * (Dot(v, P - A) / Dot(v, v));
    47	} // 点在直线上投影

    48	struct Circle {
    49	    Point c;
    50	    double r;
    51	    Circle()
    52	        : c(Point(0, 0))
    53	        , r(0) {
    54	    }
    55	    Circle(Point c, double r = 0)
    56	        : c(c)
    57	        , r(r) {
    58	    }
    59	    Point point(double a) {
    60	        return Point(c.x + cos(a) * r, c.y + sin(a) * r);
    61	    } // 输入极角返回点坐标
    62	};    // 圆

    63	// a[i] 和 b[i] 分别是第i条切线在圆A和圆B上的切点
    64	int getTangents(Circle A, Circle B, Point *a, Point *b) {
    65	    int cnt = 0;
    66	    if (A.r < B.r) {
    67	        swap(A, B);
    68	        swap(a, b);
    69	    }
    70	    double d2 =
    71	        (A.c.x - B.c.x) * (A.c.x - B.c.x) + (A.c.y - B.c.y) * (A.c.y - B.c.y);
    72	    double rdiff = A.r - B.r;
    73	    double rsum = A.r + B.r;
    74	    if (dcmp(d2 - rdiff * rdiff) < 0) return 0; // 内含

    75	    double base = atan2(B.c.y - A.c.y, B.c.x - A.c.x);
    76	    if (dcmp(d2) == 0 && dcmp(A.r - B.r) == 0) return -1; // 无限多条切线
    77	    if (dcmp(d2 - rdiff * rdiff) == 0) { // 内切，一条切线
    78	        a[cnt] = A.point(base);
    79	        b[cnt] = B.point(base);
    80	        ++cnt;
    81	        return 1;
    82	    }
    83	    // 有外公切线
    84	    double ang = acos(rdiff / sqrt(d2));
    85	    a[cnt] = A.point(base + ang);
    86	    b[cnt] = B.point(base + ang);
    87	    ++cnt;
    88	    a[cnt] = A.point(base - ang);
    89	    b[cnt] = B.point(base - ang);
    90	    ++cnt;
    91	    if (dcmp(d2 - rsum * rsum) == 0) { // 一条内公切线
    92	        a[cnt] = A.point(base);
    93	        b[cnt] = B.point(PI + base);
    94	        ++cnt;
    95	    } else if (dcmp(d2 - rsum * rsum) > 0) { // 两条内公切线
    96	        double ang = acos(rsum / sqrt(d2));
    97	        a[cnt] = A.point(base + ang);
    98	        b[cnt] = B.point(PI + base + ang);
    99	        ++cnt;
   100	        a[cnt] = A.point(base - ang);
   101	        b[cnt] = B.point(PI + base - ang);
   102	        ++cnt;
   103	    }
   104	    return cnt;
   105	} // 两圆公切线 返回切线的条数，-1表示无穷多条切线

   106	Circle Inversion_C2C(Point O, double R, Circle A) {
   107	    double OA = Length(A.c - O);
   108	    double RB = 0.5 * ((1 / (OA - A.r)) - (1 / (OA + A.r))) * R * R;
   109	    double OB = OA * RB / A.r;
   110	    double Bx = O.x + (A.c.x - O.x) * OB / OA;
   111	    double By = O.y + (A.c.y - O.y) * OB / OA;
   112	    return Circle(Point(Bx, By), RB);
   113	} // 点 O 在圆 A 外，求圆 A 的反演圆 B，R 是反演半径

   114	Circle Inversion_L2C(Point O, double R, Point A, Vector v) {
   115	    Point P = GetLineProjection(O, A, A + v);
   116	    double d = Length(O - P);
   117	    double RB = R * R / (2 * d);
   118	    Vector VB = (P - O) / d * RB;
   119	    return Circle(O + VB, RB);
   120	} // 直线反演为过 O 点的圆 B，R 是反演半径

   121	bool theSameSideOfLine(Point A, Point B, Point S, Vector v) {
   122	    return dcmp(Cross(A - S, v)) * dcmp(Cross(B - S, v)) > 0;
   123	} // 返回 true 如果 A B 两点在直线同侧

   124	int main() {
   125	    int T;
   126	    scanf("%d", &T);
   127	    while (T--) {
   128	        Circle A, B;
   129	        Point P;
   130	        scanf("%lf%lf%lf", &A.c.x, &A.c.y, &A.r);
   131	        scanf("%lf%lf%lf", &B.c.x, &B.c.y, &B.r);
   132	        scanf("%lf%lf", &P.x, &P.y);
   133	        Circle NA = Inversion_C2C(P, 10, A);
   134	        Circle NB = Inversion_C2C(P, 10, B);
   135	        Point LA[N], LB[N];
   136	        Circle ansC[N];
   137	        int q = getTangents(NA, NB, LA, LB), ans = 0;
   138	        for (int i = 0; i < q; ++i)
   139	            if (theSameSideOfLine(NA.c, NB.c, LA[i], LB[i] - LA[i])) {
   140	                if (!theSameSideOfLine(P, NA.c, LA[i], LB[i] - LA[i])) continue;
   141	                ansC[ans++] = Inversion_L2C(P, 10, LA[i], LB[i] - LA[i]);
   142	            }
   143	        printf("%d\n", ans);
   144	        for (int i = 0; i < ans; ++i) {
   145	            printf("%.8f %.8f %.8f\n", ansC[i].c.x, ansC[i].c.y, ansC[i].r);
   146	        }
   147	    }

   148	    return 0;
   149	}
 
```
 
```cpp
//./steiner_tree/luogu4294.cpp
     1	#include <bits/stdc++.h>

     2	using namespace std;

     3	#define mp make_pair
     4	typedef pair<int, int> P;
     5	typedef pair<P, int> PP;
     6	const int INF = 0x3f3f3f3f;
     7	const int dx[] = {0, 0, -1, 1};
     8	const int dy[] = {1, -1, 0, 0};
     9	int n, m, K, root;
    10	int f[101][1111], a[101], ans[11][11];
    11	bool inq[101];
    12	PP pre[101][1111];
    13	queue<P> q;

    14	bool legal(P u) {
    15	    if (u.first >= 0 && u.second >= 0 && u.first < n && u.second < m) {
    16	        return true;
    17	    }
    18	    return false;
    19	}

    20	int num(P u) {
    21	    return u.first * m + u.second;
    22	}

    23	void spfa(int s) {
    24	    memset(inq, 0, sizeof(inq));
    25	    while (!q.empty()) {
    26	        P u = q.front();
    27	        q.pop();
    28	        inq[num(u)] = 0;
    29	        for (int d = 0; d < 4; d++) {
    30	            P v = mp(u.first + dx[d], u.second + dy[d]);
    31	            int du = num(u), dv = num(v);
    32	            if (legal(v) && f[dv][s] > f[du][s] + a[dv]) {
    33	                f[dv][s] = f[du][s] + a[dv];
    34	                if (!inq[dv]) {
    35	                    inq[dv] = 1;
    36	                    q.push(v);
    37	                }
    38	                pre[dv][s] = mp(u, s);
    39	            }
    40	        }
    41	    }
    42	}

    43	void dfs(P u, int s) {
    44	    if (!pre[num(u)][s].second) return;
    45	    ans[u.first][u.second] = 1;
    46	    int nu = num(u);
    47	    if (pre[nu][s].first == u) dfs(u, s ^ pre[nu][s].second);
    48	    dfs(pre[nu][s].first, pre[nu][s].second);
    49	}

    50	int main() {
    51	    memset(f, INF, sizeof(f));
    52	    scanf("%d %d", &n, &m);
    53	    int tot = 0;
    54	    for (int i = 0; i < n; i++) {
    55	        for (int j = 0; j < m; j++) {
    56	            scanf("%d", &a[tot]);
    57	            if (!a[tot]) {
    58	                f[tot][1 << (K++)] = 0;
    59	                root = tot;
    60	            }
    61	            tot++;
    62	        }
    63	    }
    64	    for (int s = 1; s < (1 << K); s++) {
    65	        for (int i = 0; i < n * m; i++) {
    66	            for (int subs = s & (s - 1); subs; subs = s & (subs - 1)) {
    67	                if (f[i][s] > f[i][subs] + f[i][s ^ subs] - a[i]) {
    68	                    f[i][s] = f[i][subs] + f[i][s ^ subs] - a[i];
    69	                    pre[i][s] = mp(mp(i / m, i % m), subs);
    70	                }
    71	            }
    72	            if (f[i][s] < INF) q.push(mp(i / m, i % m));
    73	        }
    74	        spfa(s);
    75	    }
    76	    printf("%d\n", f[root][(1 << K) - 1]);
    77	    dfs(mp(root / m, root % m), (1 << K) - 1);
    78	    for (int i = 0, tot = 0; i < n; i++) {
    79	        for (int j = 0; j < m; j++) {
    80	            if (!a[tot++])
    81	                putchar('x');
    82	            else
    83	                putchar(ans[i][j] ? 'o' : '_');
    84	        }
    85	        if (i != n - 1) printf("\n");
    86	    }
    87	    return 0;
    88	}
 
```
 
```cpp
//./steiner_tree/luogu6192.cpp
     1	#include <bits/stdc++.h>

     2	using namespace std;

     3	const int maxn = 510;
     4	const int INF = 0x3f3f3f3f;
     5	typedef pair<int, int> P;
     6	int n, m, k;

     7	struct edge {
     8	    int to, next, w;
     9	} e[maxn << 1];

    10	int head[maxn << 1], tree[maxn << 1], tot;
    11	int dp[maxn][5000], vis[maxn];
    12	int key[maxn];
    13	priority_queue<P, vector<P>, greater<P>> q;

    14	void add(int u, int v, int w) {
    15	    e[++tot] = edge{v, head[u], w};
    16	    head[u] = tot;
    17	}

    18	void dijkstra(int s) {
    19	    memset(vis, 0, sizeof(vis));
    20	    while (!q.empty()) {
    21	        P item = q.top();
    22	        q.pop();
    23	        if (vis[item.second]) continue;
    24	        vis[item.second] = 1;
    25	        for (int i = head[item.second]; i; i = e[i].next) {
    26	            if (dp[tree[i]][s] > dp[item.second][s] + e[i].w) {
    27	                dp[tree[i]][s] = dp[item.second][s] + e[i].w;
    28	                q.push(P(dp[tree[i]][s], tree[i]));
    29	            }
    30	        }
    31	    }
    32	}

    33	/*
    34	给定一个包含 n 个结点和 m 条带权边的无向连通图 G=(V,E)。

    35	再给定包含 k 个结点的点集 S，选出 G 的子图
    36	G'=(V',E')，使得：
    37	    S⊆V′；
    38	    G′ 为连通图；
    39	    E′ 中所有边的权值和最小。
    40	你只需要求出 E′ 中所有边的权值和。
    41	*/
    42	int main() {
    43	    memset(dp, INF, sizeof(dp));
    44	    scanf("%d %d %d", &n, &m, &k);
    45	    int u, v, w;
    46	    for (int i = 1; i <= m; i++) {
    47	        scanf("%d %d %d", &u, &v, &w);
    48	        add(u, v, w);
    49	        tree[tot] = v;
    50	        add(v, u, w);
    51	        tree[tot] = u;
    52	    }
    53	    for (int i = 1; i <= k; i++) {
    54	        scanf("%d", &key[i]);
    55	        dp[key[i]][1 << (i - 1)] = 0;
    56	    }
    57	    for (int s = 1; s < (1 << k); s++) {
    58	        for (int i = 1; i <= n; i++) {
    59	            for (int subs = s & (s - 1); subs; subs = s & (subs - 1))
    60	                dp[i][s] = min(dp[i][s], dp[i][subs] + dp[i][s ^ subs]);
    61	            if (dp[i][s] != INF) q.push(P(dp[i][s], i));
    62	        }
    63	        dijkstra(s);
    64	    }
    65	    printf("%d\n", dp[key[1]][(1 << k) - 1]);
    66	    return 0;
    67	}

 
```
 
```cpp
//./dbg_func/dbg_min.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;

     3	#define show(arg) cerr << "\033[36m" << #arg << "\033[0m = " << arg << endl;
     4	#define dbg(...) debug(string(#__VA_ARGS__).begin(), __VA_ARGS__)

     5	void debug(string::iterator it) {
     6	    cerr << endl;
     7	}

     8	void debug(string::iterator it, auto a, auto... args) {
     9	    cerr << "\033[36m";
    10	    int cnt = 0;
    11	    while (*it != ',' || cnt) {
    12	        if (*it == '(' || *it == '{') cnt++;
    13	        if (*it == ')' || *it == '}') cnt--;
    14	        if (*it != ' ') cerr << *it;
    15	        ++it;
    16	    }
    17	    cerr << "\033[0m"
    18	         << " = " << a << "   ";
    19	    debug(++it, args...);
    20	}
 
```
 
```cpp
//./intersection_of_half_planes/intersection_of_half_planes.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	#define ll long long
     4	const double eps = 1e-6;
     5	const int maxn = 1e5 + 7;

     6	// tips: this tmp maybe need each vector has intersection

     7	bool dcmp(double a, double b) {
     8	    return fabs(a - b) < eps;
     9	}

    10	class point
    11	// 点或向量的类
    12	{
    13	  public:
    14	    double x, y;
    15	    void read() {
    16	        cin >> x >> y;
    17	    }
    18	    void print() {
    19	        cout << x << y;
    20	    }
    21	    point operator+(const point &p) const {
    22	        return {x + p.x, y + p.y};
    23	    }
    24	    point operator-(const point &p) const {
    25	        return {x - p.x, y - p.y};
    26	    }
    27	    point operator*(double p) const {
    28	        return {x * p, y * p};
    29	    }
    30	    point operator/(double p) const {
    31	        return {x / p, y / p};
    32	    }
    33	    double operator*(const point &p) const {
    34	        return x * p.x + y * p.y;
    35	    }
    36	    double len() {
    37	        return sqrt(x * x + y * y);
    38	    }
    39	    double dis(const point &p) const {
    40	        return ((*this) - p).len();
    41	    }
    42	    double angle() {
    43	        return atan2(y, x);
    44	    }
    45	} ans[maxn], tmpp[maxn];

    46	double cross(point a, point b)
    47	//求向量外积
    48	{
    49	    return a.x * b.y - a.y * b.x;
    50	}

    51	double area(point s[], int n)
    52	//求凸多边形面积
    53	{
    54	    double ret = 0;
    55	    s[n + 1] = s[1];
    56	    for (int i = 1; i <= n; i++)
    57	        ret += cross(s[i], s[i + 1]);
    58	    return fabs(ret / 2);
    59	}

    60	struct line
    61	// 线
    62	{
    63	    point s, t;
    64	    double ang;
    65	    void getline(point a, point b) {
    66	        s = a, t = b;
    67	        ang = (t - s).angle();
    68	    }
    69	    line() {
    70	    }
    71	    line(point a, point b) {
    72	        this->s = a, this->t = b, this->ang = (a - b).angle();
    73	    }
    74	} li[maxn], tmpl[maxn];

    75	point intersection(line a, line b)
    76	// 交点
    77	{
    78	    return a.s + (a.t - a.s) * (cross(a.s - b.s, b.t - b.s) /
    79	                                cross(b.t - b.s, a.t - a.s));
    80	}

    81	bool parallel(line a, line b)
    82	//平行
    83	{
    84	    return dcmp(cross(a.t - a.s, b.t - b.s), 0);
    85	}

    86	bool isrig(line a, point b)
    87	//判断点是否在向量右边
    88	{
    89	    return cross(a.t - a.s, b - a.s) + eps < 0;
    90	}

    91	bool linecmp(line a, line b)
    92	// 极角序
    93	{
    94	    return dcmp(a.ang, b.ang) ? cross(a.t - a.s, b.t - a.s) + eps < 0
    95	                              : a.ang < b.ang;
    96	}

    97	bool SI(line *li, int n, point *ret, int &m, line *ql, point *qp)
    98	// 半平面交，li是向量集合，n是数量，ret是答案，ql和qp是临时的
    99	{
   100	    int l, r;
   101	    sort(li + 1, li + 1 + n, linecmp);
   102	    ql[l = r = 1] = li[1];
   103	    for (int i = 2; i <= n; i++)
   104	        if (!dcmp(li[i].ang, li[i - 1].ang)) {
   105	            while (l < r && isrig(li[i], qp[r - 1]))
   106	                --r;
   107	            while (l < r && isrig(li[i], qp[l]))
   108	                ++l;
   109	            if (l <= r) qp[r] = intersection(ql[r], li[i]);
   110	            ql[++r] = li[i];
   111	            if (l < r &&
   112	                (parallel(ql[l], ql[l + 1]) || parallel(ql[r], ql[r - 1])))
   113	                return false;
   114	        }
   115	    while (l < r && isrig(ql[l], qp[r - 1]))
   116	        --r;
   117	    while (l < r && isrig(ql[r], qp[l]))
   118	        ++l;
   119	    if (r - l <= 1) return false;
   120	    qp[r] = intersection(ql[l], ql[r]);
   121	    m = 0;
   122	    for (int i = l; i <= r; i++)
   123	        ret[++m] = qp[i];
   124	    return true;
   125	}

   126	point a[maxn], b[maxn];

   127	int main() {
   128	    int t;
   129	    cin >> t;
   130	    int tot = 0;
   131	    while (t--) {
   132	        int n;
   133	        cin >> n;
   134	        for (int i = 0; i < n; i++)
   135	            a[i].read();
   136	        for (int i = 0; i < n; i++)
   137	            li[++tot].getline(a[i], a[(i + 1) % n]);
   138	    }
   139	    int m = 0;
   140	    if (SI(li, tot, ans, m, tmpl, tmpp))
   141	        printf("%.3f", area(ans, m));
   142	    else
   143	        cout << "0.000";
   144	    return 0;
   145	}
 
```
 
```cpp
//./LGV/hdu5852.cpp
     1	#include <bits/stdc++.h>
     2	typedef long long ll;
     3	const int K = 105;
     4	const int N = 100005;
     5	const int mod = 1e9 + 7;
     6	int T, n, k, a[K], b[K], fact[N << 1], m[K][K];
     7	int qpow(int x, int y) {
     8	    int out = 1;
     9	    while (y) {
    10	        if (y & 1) out = (ll) out * x % mod;
    11	        x = (ll) x * x % mod;
    12	        y >>= 1;
    13	    }
    14	    return out;
    15	}
    16	int c(int x, int y) {
    17	    return (ll) fact[x] * qpow(fact[y], mod - 2) % mod *
    18	           qpow(fact[x - y], mod - 2) % mod;
    19	}
    20	// 有一个 nxn 的棋盘，一个棋子从 $(x, y)$ 只能走到 $(x, y+1)$ 或 $(x +
    21	// 1, y)$ ，有 $k$ 个棋子，一开始第 $i$ 个棋子放在 $(1, a_i)$ ，最终要到 $(n,
    22	// b_i)$ ，路径要两两不相交，求方案数对 $10^9+7$ 取模。
    23	int main() {
    24	    fact[0] = 1;
    25	    for (int i = 1; i < N * 2; ++i)
    26	        fact[i] = (ll) fact[i - 1] * i % mod;
    27	    scanf("%d", &T);
    28	    while (T--) {
    29	        scanf("%d%d", &n, &k);
    30	        for (int i = 1; i <= k; ++i)
    31	            scanf("%d", a + i);
    32	        for (int i = 1; i <= k; ++i)
    33	            scanf("%d", b + i);
    34	        for (int i = 1; i <= k; ++i) {
    35	            for (int j = 1; j <= k; ++j) {
    36	                if (a[i] <= b[j])
    37	                    m[i][j] = c(b[j] - a[i] + n - 1, n - 1);
    38	                else
    39	                    m[i][j] = 0;
    40	            }
    41	        }
    42	        for (int i = 1; i < k; ++i) {
    43	            if (!m[i][i]) {
    44	                for (int j = i + 1; j <= k; ++j) {
    45	                    if (m[j][i]) {
    46	                        std::swap(m[i], m[j]);
    47	                        break;
    48	                    }
    49	                }
    50	            }
    51	            if (!m[i][i]) continue;
    52	            int inv = qpow(m[i][i], mod - 2);
    53	            for (int j = i + 1; j <= k; ++j) {
    54	                if (!m[j][i]) continue;
    55	                int mul = (ll) m[j][i] * inv % mod;
    56	                for (int p = i; p <= k; ++p) {
    57	                    m[j][p] = (m[j][p] - (ll) m[i][p] * mul % mod + mod) % mod;
    58	                }
    59	            }
    60	        }
    61	        int ans = 1;
    62	        for (int i = 1; i <= k; ++i)
    63	            ans = (ll) ans * m[i][i] % mod;
    64	        printf("%d\n", ans);
    65	    }
    66	    return 0;
    67	}
 
```
 
```cpp
//./palindromic_tree/luogu3649.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	typedef long long ll;
     4	const int maxn = 300000 + 5;
     5	namespace pam {
     6	int sz, tot, last;
     7	int cnt[maxn], ch[maxn][26], len[maxn], fail[maxn];
     8	char s[maxn];
     9	int node(int l) {
    10	    sz++;
    11	    memset(ch[sz], 0, sizeof(ch[sz]));
    12	    len[sz] = l;
    13	    fail[sz] = cnt[sz] = 0;
    14	    return sz;
    15	}
    16	void clear() {
    17	    sz = -1;
    18	    last = 0;
    19	    s[tot = 0] = '$';
    20	    node(0);
    21	    node(-1);
    22	    fail[0] = 1;
    23	}
    24	int getfail(int x) {
    25	    while (s[tot - len[x] - 1] != s[tot])
    26	        x = fail[x];
    27	    return x;
    28	}
    29	void insert(char c) {
    30	    s[++tot] = c;
    31	    int now = getfail(last);
    32	    if (!ch[now][c - 'a']) {
    33	        int x = node(len[now] + 2);
    34	        fail[x] = ch[getfail(fail[now])][c - 'a'];
    35	        ch[now][c - 'a'] = x;
    36	    }
    37	    last = ch[now][c - 'a'];
    38	    cnt[last]++;
    39	}
    40	ll solve() {
    41	    ll ans = 0;
    42	    for (int i = sz; i >= 0; i--) {
    43	        cnt[fail[i]] += cnt[i];
    44	    }
    45	    for (int i = 1; i <= sz; i++) {
    46	        ans = max(ans, 1ll * len[i] * cnt[i]);
    47	    }
    48	    return ans;
    49	}
    50	} // namespace pam
    51	char s[maxn];
    52	// 定义 $s$ 的一个子串的存在值为这个子串在 $s$
    53	// 中出现的次数乘以这个子串的长度。对于给定的字符串 $s$
    54	// ，求所有回文子串中的最大存在值。
    55	int main() {
    56	    pam::clear();
    57	    scanf("%s", s + 1);
    58	    for (int i = 1; s[i]; i++) {
    59	        pam::insert(s[i]);
    60	    }
    61	    printf("%lld\n", pam::solve());
    62	    return 0;
    63	}

 
```
 
```cpp
//./palindromic_tree/CF932G.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	typedef long long ll;
     4	const int mod = 1e9 + 7;
     5	const int maxn = 1000000 + 5;
     6	inline int add(int x, int y) {
     7	    x += y;
     8	    return x >= mod ? x -= mod : x;
     9	}
    10	namespace pam {
    11	int sz, tot, last;
    12	int ch[maxn][26], len[maxn], fail[maxn];
    13	int cnt[maxn], dep[maxn], dif[maxn], slink[maxn];
    14	char s[maxn];
    15	int node(int l) {
    16	    sz++;
    17	    memset(ch[sz], 0, sizeof(ch[sz]));
    18	    len[sz] = l;
    19	    fail[sz] = 0;
    20	    cnt[sz] = 0;
    21	    dep[sz] = 0;
    22	    return sz;
    23	}
    24	void clear() {
    25	    sz = -1;
    26	    last = 0;
    27	    s[tot = 0] = '$';
    28	    node(0);
    29	    node(-1);
    30	    fail[0] = 1;
    31	}
    32	int getfail(int x) {
    33	    while (s[tot - len[x] - 1] != s[tot])
    34	        x = fail[x];
    35	    return x;
    36	}
    37	void insert(char c) {
    38	    s[++tot] = c;
    39	    int now = getfail(last);
    40	    if (!ch[now][c - 'a']) {
    41	        int x = node(len[now] + 2);
    42	        fail[x] = ch[getfail(fail[now])][c - 'a'];
    43	        dep[x] = dep[fail[x]] + 1;
    44	        ch[now][c - 'a'] = x;
    45	        dif[x] = len[x] - len[fail[x]];
    46	        if (dif[x] == dif[fail[x]])
    47	            slink[x] = slink[fail[x]];
    48	        else
    49	            slink[x] = fail[x];
    50	    }
    51	    last = ch[now][c - 'a'];
    52	    cnt[last]++;
    53	}
    54	} // namespace pam
    55	using pam::dif;
    56	using pam::fail;
    57	using pam::len;
    58	using pam::slink;
    59	int n, dp[maxn], g[maxn];
    60	char s[maxn], t[maxn];
    61	int main() {
    62	    pam::clear();
    63	    scanf("%s", s + 1);
    64	    n = strlen(s + 1);
    65	    for (int i = 1, j = 0; i <= n; i++)
    66	        t[++j] = s[i], t[++j] = s[n - i + 1];
    67	    dp[0] = 1;
    68	    for (int i = 1; i <= n; i++) {
    69	        pam::insert(t[i]);
    70	        for (int x = pam::last; x > 1; x = slink[x]) {
    71	            g[x] = dp[i - len[slink[x]] - dif[x]];
    72	            if (dif[x] == dif[fail[x]]) g[x] = add(g[x], g[fail[x]]);
    73	            if (i % 2 == 0) dp[i] = add(dp[i], g[x]);
    74	        }
    75	    }
    76	    printf("%d", dp[n]);
    77	    return 0;
    78	}

 
```
 
```cpp
//./segment_tree_beats/hdu5306.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	#define N 500010
     4	#define LL long long
     5	int T, n, m, a[N];
     6	#define lc (p << 1)
     7	#define rc (p << 1 | 1)
     8	struct tree {
     9	    int l, r, mx, sx, cx;
    10	    LL sum;
    11	} t[N << 2];
    12	void up(int p) {
    13	    t[p].cx = 0;
    14	    t[p].sum = t[lc].sum + t[rc].sum;
    15	    t[p].mx = max(t[lc].mx, t[rc].mx);
    16	    t[p].sx = max(t[lc].sx, t[rc].sx);
    17	    if (t[lc].mx ^ t[rc].mx) t[p].sx = max(t[p].sx, min(t[lc].mx, t[rc].mx));
    18	    if (t[lc].mx == t[p].mx) t[p].cx += t[lc].cx;
    19	    if (t[rc].mx == t[p].mx) t[p].cx += t[rc].cx;
    20	}
    21	void build(int p, int l, int r) {
    22	    t[p].l = l;
    23	    t[p].r = r;
    24	    if (l == r) {
    25	        t[p].sum = t[p].mx = a[l];
    26	        t[p].sx = -1;
    27	        t[p].cx = 1;
    28	        return;
    29	    }
    30	    int mid = (l + r) >> 1;
    31	    build(lc, l, mid);
    32	    build(rc, mid + 1, r);
    33	    up(p);
    34	}
    35	void MIN(int p, int v) {
    36	    if (v >= t[p].mx) return;
    37	    t[p].sum += 1ll * (v - t[p].mx) * t[p].cx;
    38	    t[p].mx = v;
    39	}
    40	void down(int p) {
    41	    MIN(lc, t[p].mx);
    42	    MIN(rc, t[p].mx);
    43	}
    44	void Min(int p, int l, int r, int v) {
    45	    if (v >= t[p].mx) return;
    46	    if (l <= t[p].l && t[p].r <= r && t[p].sx < v) {
    47	        MIN(p, v);
    48	        return;
    49	    }
    50	    int mid = (t[p].l + t[p].r) >> 1;
    51	    down(p);
    52	    if (l <= mid) Min(lc, l, r, v);
    53	    if (r > mid) Min(rc, l, r, v);
    54	    up(p);
    55	}
    56	LL query(int p, int l, int r, int op) {
    57	    if (l <= t[p].l && t[p].r <= r) {
    58	        if (op == 2) return t[p].sum;
    59	        return t[p].mx;
    60	    }
    61	    int mid = (t[p].l + t[p].r) >> 1;
    62	    down(p);
    63	    LL ans;
    64	    if (op == 2) {
    65	        ans = 0;
    66	        if (l <= mid) ans += query(lc, l, r, op);
    67	        if (r > mid) ans += query(rc, l, r, op);
    68	    } else {
    69	        ans = -1;
    70	        if (l <= mid) ans = max(ans, query(lc, l, r, op));
    71	        if (r > mid) ans = max(ans, query(rc, l, r, op));
    72	    }
    73	    up(p);
    74	    return ans;
    75	}
    76	// 给一个序列，要求支持区间取min（即对于一段区间，用min(a[i],x)替换a[i]（x已给出）），询问区间和以及区间最大值。
    77	// command：0 替换 1 求最大值 2 求和
    78	int main() {
    79	    int op, x, y, z;
    80	    scanf("%d", &T);
    81	    while (T--) {
    82	        scanf("%d%d", &n, &m);
    83	        for (int i = 1; i <= n; i++)
    84	            scanf("%d", &a[i]);
    85	        build(1, 1, n);
    86	        for (int i = 1; i <= m; i++) {
    87	            scanf("%d%d%d", &op, &x, &y);
    88	            if (!op) {
    89	                scanf("%d", &z);
    90	                Min(1, x, y, z);
    91	            } else {
    92	                printf("%lld\n", query(1, x, y, op));
    93	            }
    94	        }
    95	    }
    96	    return 0;
    97	}
 
```
 
```cpp
//./segment_tree_beats/bzoj4965.cpp
     1	#include <bits/stdc++.h>
     2	using namespace std;
     3	#define N 500010
     4	#define LL long long
     5	#define INF 1 << 30
     6	int n, m, a[N];
     7	#define lc (p << 1)
     8	#define rc (p << 1 | 1)
     9	struct node {
    10	    int l, r, len, mx, mn, sx, sn, cx, cn;
    11	    LL laz, sum;
    12	} t[N << 2];
    13	void up(int p) {
    14	    t[p].cx = t[p].cn = 0;
    15	    t[p].sum = t[lc].sum + t[rc].sum;
    16	    t[p].mx = max(t[lc].mx, t[rc].mx);
    17	    t[p].sx = max(t[lc].sx, t[rc].sx);
    18	    if (t[lc].mx ^ t[rc].mx) t[p].sx = max(t[p].sx, min(t[lc].mx, t[rc].mx));
    19	    if (t[p].mx == t[lc].mx) t[p].cx += t[lc].cx;
    20	    if (t[p].mx == t[rc].mx) t[p].cx += t[rc].cx;
    21	    t[p].mn = min(t[lc].mn, t[rc].mn);
    22	    t[p].sn = min(t[lc].sn, t[rc].sn);
    23	    if (t[lc].mn ^ t[rc].mn) t[p].sn = min(t[p].sn, max(t[lc].mn, t[rc].mn));
    24	    if (t[p].mn == t[lc].mn) t[p].cn += t[lc].cn;
    25	    if (t[p].mn == t[rc].mn) t[p].cn += t[rc].cn;
    26	}
    27	void build(int p, int l, int r) {
    28	    t[p].l = l;
    29	    t[p].r = r;
    30	    t[p].len = r - l + 1;
    31	    if (l == r) {
    32	        t[p].mx = t[p].mn = t[p].sum = a[l];
    33	        t[p].cx = t[p].cn = 1;
    34	        t[p].sx = -INF;
    35	        t[p].sn = INF;
    36	        return;
    37	    }
    38	    int mid = (l + r) >> 1;
    39	    build(lc, l, mid);
    40	    build(rc, mid + 1, r);
    41	    up(p);
    42	}
    43	void now(int p, int v) {
    44	    t[p].sum += 1ll * t[p].len * v;
    45	    t[p].laz += v;
    46	    t[p].mx += v;
    47	    t[p].mn += v;
    48	    if (t[p].sx != -INF) t[p].sx += v;
    49	    if (t[p].sn != INF) t[p].sn += v;
    50	}
    51	void MAX(int p, int v) {
    52	    t[p].sum += 1ll * (v - t[p].mn) * t[p].cn;
    53	    t[p].mn = v;
    54	    t[p].mx = max(v, t[p].mx);
    55	    if (t[p].mx == t[p].mn) {
    56	        t[p].sum = 1ll * t[p].len * v;
    57	        t[p].cn = t[p].cx = t[p].len;
    58	        t[p].sx = -INF;
    59	        t[p].sn = INF;
    60	    } else
    61	        t[p].sx = max(v, t[p].sx);
    62	}
    63	void MIN(int p, int v) {
    64	    t[p].sum += 1ll * (v - t[p].mx) * t[p].cx;
    65	    t[p].mx = v;
    66	    t[p].mn = min(v, t[p].mn);
    67	    if (t[p].mx == t[p].mn) {
    68	        t[p].sum = 1ll * t[p].len * v;
    69	        t[p].cn = t[p].cx = t[p].len;
    70	        t[p].sx = -INF;
    71	        t[p].sn = INF;
    72	    } else
    73	        t[p].sn = min(v, t[p].sn);
    74	}
    75	void down(int p) {
    76	    if (t[p].laz) {
    77	        now(lc, t[p].laz);
    78	        now(rc, t[p].laz);
    79	        t[p].laz = 0;
    80	    }
    81	    if (t[lc].mn < t[p].mn && t[p].mn < t[lc].sn) MAX(lc, t[p].mn);
    82	    if (t[rc].mn < t[p].mn && t[p].mn < t[rc].sn) MAX(rc, t[p].mn);
    83	    if (t[lc].sx < t[p].mx && t[p].mx < t[lc].mx) MIN(lc, t[p].mx);
    84	    if (t[rc].sx < t[p].mx && t[p].mx < t[rc].mx) MIN(rc, t[p].mx);
    85	}
    86	void add(int p, int l, int r, int v) {
    87	    if (!v) return;
    88	    if (l <= t[p].l && t[p].r <= r) {
    89	        now(p, v);
    90	        return;
    91	    }
    92	    int mid = (t[p].l + t[p].r) >> 1;
    93	    down(p);
    94	    if (l <= mid) add(lc, l, r, v);
    95	    if (r > mid) add(rc, l, r, v);
    96	    up(p);
    97	}
    98	void Max(int p, int l, int r, int v) {
    99	    if (t[p].mn >= v) return;
   100	    if (l <= t[p].l && t[p].r <= r && v < t[p].sn) {
   101	        MAX(p, v);
   102	        return;
   103	    }
   104	    int mid = (t[p].l + t[p].r) >> 1;
   105	    down(p);
   106	    if (l <= mid) Max(lc, l, r, v);
   107	    if (r > mid) Max(rc, l, r, v);
   108	    up(p);
   109	}
   110	void Min(int p, int l, int r, int v) {
   111	    if (t[p].mx <= v) return;
   112	    if (l <= t[p].l && t[p].r <= r && v > t[p].sx) {
   113	        MIN(p, v);
   114	        return;
   115	    }
   116	    int mid = (t[p].l + t[p].r) >> 1;
   117	    down(p);
   118	    if (l <= mid) Min(lc, l, r, v);
   119	    if (r > mid) Min(rc, l, r, v);
   120	    up(p);
   121	}
   122	LL query(int p, int l, int r, int op) {
   123	    if (l <= t[p].l && t[p].r <= r) {
   124	        if (op == 4) return t[p].sum;
   125	        if (op == 5) return t[p].mx;
   126	        if (op == 6) return t[p].mn;
   127	    }
   128	    int mid = (t[p].l + t[p].r) >> 1;
   129	    LL ans;
   130	    down(p);
   131	    if (op == 4) {
   132	        ans = 0;
   133	        if (l <= mid) ans += query(lc, l, r, op);
   134	        if (r > mid) ans += query(rc, l, r, op);
   135	    }
   136	    if (op == 5) {
   137	        ans = -INF;
   138	        if (l <= mid) ans = max(ans, query(lc, l, r, op));
   139	        if (r > mid) ans = max(ans, query(rc, l, r, op));
   140	    }
   141	    if (op == 6) {
   142	        ans = INF;
   143	        if (l <= mid) ans = min(ans, query(lc, l, r, op));
   144	        if (r > mid) ans = min(ans, query(rc, l, r, op));
   145	    }
   146	    up(p);
   147	    return ans;
   148	}
   149	/*
   150	    1.给一个区间 [L, R] 加上一个数x
   151	　　2.把一个区间 [L, R] 里小于 x 的数变成x
   152	　　3.把一个区间 [L, R] 里大于 x 的数变成x
   153	　　4.求区间 [L, R] 的和
   154	　　5.求区间 [L, R] 的最大值
   155	　　6.求区间 [L, R] 的最小值
   156	*/
   157	int main() {
   158	    int op, x, y, z;
   159	    scanf("%d", &n);
   160	    for (int i = 1; i <= n; i++)
   161	        scanf("%d", &a[i]);
   162	    scanf("%d", &m);
   163	    build(1, 1, n);
   164	    for (int i = 1; i <= m; i++) {
   165	        scanf("%d%d%d", &op, &x, &y);
   166	        if (op <= 3) {
   167	            scanf("%d", &z);
   168	            if (op == 1) add(1, x, y, z);
   169	            if (op == 2) Max(1, x, y, z);
   170	            if (op == 3) Min(1, x, y, z);
   171	        } else {
   172	            printf("%lld\n", query(1, x, y, op));
   173	        }
   174	    }
   175	    return 0;
   176	}
 
```
