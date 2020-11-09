#include <bits/extc++.h>
#include <bits/stdc++.h>
/*
若编译器没有 extc++.h 使用如下三个库
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/priority_queue.hpp>
#include <ext/pb_ds/tree_policy.hpp>
*/
class T {

  public:
    int val, rnk;
    T(int val, int rnk) {
        this->val = val;
        this->rnk = rnk;
    }
    bool operator<(const T &b) const { return this->val < b.val; }
    bool operator>(const T &b) const {
        if (this->val != b.val)
            return this->val > b.val;
        return this->rnk < b.rnk;
    }
};
#define param T, greater<T>, pairing_heap_tag
using namespace std;
using namespace __gnu_pbds;
const int maxn = 1e5 + 7;
__gnu_pbds::priority_queue<param> qs[maxn];
int a[maxn];

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

int main() {
    int n, m;
    cin >> n >> m;
    s.init(n);

    for (int i = 1; i <= n; i++) {
        int val;
        cin >> val;
        qs[i].push(T(val, i));
        a[i] = val;
    }
    for (int i = 1; i <= m; i++) {
        int cmd;
        cin >> cmd;
        if (cmd == 1) {
            int x, y;
            cin >> x >> y;
            if (a[x] == -1 || a[y] == -1)
                continue;
            if (x == y)
                continue;
            qs[s.find(x)].join(qs[s.find(y)]);
            s.Union(x, y);
        } else {
            int x;
            cin >> x;
            if (a[x] == -1)
                puts("-1");
            else {
                int y = s.find(x);
                T tmp = qs[y].top();
                cout << tmp.val << endl;
                a[tmp.rnk] = -1;
                qs[y].pop();
            }
        }
    }
    return 0;
}