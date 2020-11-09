#include <bits/extc++.h>
#include <bits/stdc++.h>
#ifndef ONLINE_JUDGE
#include "dbg_func"
#else
#define dbg(...) (__VA_ARGS__)
#endif
/*
若编译器没有 extc++.h 使用如下三个库
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/priority_queue.hpp>
#include <ext/pb_ds/tree_policy.hpp>
*/
#define param int, less<int>, pairing_heap_tag
using namespace std;
using namespace __gnu_pbds;
const int maxn = 1e5 + 7;
#define int long long
/*
__gnu_pbds::priority_queue<int, greater<int>> que1;                       // （默认为配对堆）
__gnu_pbds::priority_queue<int, less<int>> que1;                          // （默认为配对堆）
__gnu_pbds::priority_queue<int, greater<int>, pairing_heap_tag> que2;     // 配对堆（应该是最快的）
__gnu_pbds::priority_queue<int, greater<int>, binary_heap_tag> que3;      // 二叉堆
__gnu_pbds::priority_queue<int, greater<int>, binomial_heap_tag> que4;    // 二项堆
__gnu_pbds::priority_queue<int, greater<int>, rc_binomial_heap_tag> que5; //
__gnu_pbds::priority_queue<int, greater<int>, thin_heap_tag> que6;        // 斐波那契堆
只有 push pop join 用 binary_heap
有 modify 用 pairing_heap & thin_heap <- thin_heap 可能会MLE，慎用
优化 Dijkstra 用 pairing_heap
std 垃圾
less 大根堆 greater 小根堆
例题 猴子选国王 并查集 + 左偏树
*/
__gnu_pbds::priority_queue<param> pq[maxn];

struct DisjointSet {
    int p[maxn];
    void init(int n) {
        for (int i = 1; i <= n; i++) {
            p[i] = i;
        }
    }
    int find(int u) { return u == p[u] ? u : p[u] = find(p[u]); }
    bool isJoined(int u, int v) { return find(u) == find(v); }

    void Union(int u, int v) {
        u = find(u);
        v = find(v);
        p[v] = u;
    }
} s;
int c[maxn], l[maxn], deg[maxn], totc[maxn];

signed main() {
    /*
    修改和删除示例
    __gnu_pbds::priority_queue<int> testplus;
    auto it = testplus.push(0);
    testplus.push(1);
    testplus.push(2);
    testplus.modify(it, 3);
    cout << testplus.top();
    testplus.erase(it);
    cout << testplus.top();
    */
    int n, m;
    cin >> n >> m;
    s.init(n);
    for (int i = 1; i <= n; i++) {
        int bi;
        cin >> bi >> c[i] >> l[i];
        s.p[i] = bi;
        deg[bi]++;
        totc[i] = c[i];
        pq[i].push(c[i]);
    }
    queue<int> ninjya;
    for (int i = 1; i <= n; i++)
        if (deg[i] == 0) {
            ninjya.push(i);
        }
    int ans = 0;
    while (ninjya.size()) {
        int now = ninjya.front();
        ninjya.pop();
        while (totc[now] > m) {
            totc[now] -= pq[now].top();
            pq[now].pop();
        }
        // cout << now << ' ' << totc[now] << endl;
        dbg(now, totc[now], l[now], pq[now].size());
        ans = max(ans, l[now] * (long long)pq[now].size());
        pq[s.p[now]].join(pq[now]);
        totc[s.p[now]] += totc[now];
        deg[s.p[now]]--;
        if (deg[s.p[now]] == 0)
            ninjya.push(s.p[now]);
    }
    cout << ans << endl;
    return 0;
}