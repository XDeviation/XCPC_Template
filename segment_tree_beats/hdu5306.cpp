#include <bits/stdc++.h>
using namespace std;
#define N 500010
#define LL long long
int T, n, m, a[N];
#define lc (p << 1)
#define rc (p << 1 | 1)
struct tree {
    int l, r, mx, sx, cx;
    LL sum;
} t[N << 2];
void up(int p) {
    t[p].cx = 0;
    t[p].sum = t[lc].sum + t[rc].sum;
    t[p].mx = max(t[lc].mx, t[rc].mx);
    t[p].sx = max(t[lc].sx, t[rc].sx);
    if (t[lc].mx ^ t[rc].mx) t[p].sx = max(t[p].sx, min(t[lc].mx, t[rc].mx));
    if (t[lc].mx == t[p].mx) t[p].cx += t[lc].cx;
    if (t[rc].mx == t[p].mx) t[p].cx += t[rc].cx;
}
void build(int p, int l, int r) {
    t[p].l = l;
    t[p].r = r;
    if (l == r) {
        t[p].sum = t[p].mx = a[l];
        t[p].sx = -1;
        t[p].cx = 1;
        return;
    }
    int mid = (l + r) >> 1;
    build(lc, l, mid);
    build(rc, mid + 1, r);
    up(p);
}
void MIN(int p, int v) {
    if (v >= t[p].mx) return;
    t[p].sum += 1ll * (v - t[p].mx) * t[p].cx;
    t[p].mx = v;
}
void down(int p) {
    MIN(lc, t[p].mx);
    MIN(rc, t[p].mx);
}
void Min(int p, int l, int r, int v) {
    if (v >= t[p].mx) return;
    if (l <= t[p].l && t[p].r <= r && t[p].sx < v) {
        MIN(p, v);
        return;
    }
    int mid = (t[p].l + t[p].r) >> 1;
    down(p);
    if (l <= mid) Min(lc, l, r, v);
    if (r > mid) Min(rc, l, r, v);
    up(p);
}
LL query(int p, int l, int r, int op) {
    if (l <= t[p].l && t[p].r <= r) {
        if (op == 2) return t[p].sum;
        return t[p].mx;
    }
    int mid = (t[p].l + t[p].r) >> 1;
    down(p);
    LL ans;
    if (op == 2) {
        ans = 0;
        if (l <= mid) ans += query(lc, l, r, op);
        if (r > mid) ans += query(rc, l, r, op);
    } else {
        ans = -1;
        if (l <= mid) ans = max(ans, query(lc, l, r, op));
        if (r > mid) ans = max(ans, query(rc, l, r, op));
    }
    up(p);
    return ans;
}
// 给一个序列，要求支持区间取min（即对于一段区间，用min(a[i],x)替换a[i]（x已给出）），询问区间和以及区间最大值。
// command：0 替换 1 求最大值 2 求和
int main() {
    int op, x, y, z;
    scanf("%d", &T);
    while (T--) {
        scanf("%d%d", &n, &m);
        for (int i = 1; i <= n; i++)
            scanf("%d", &a[i]);
        build(1, 1, n);
        for (int i = 1; i <= m; i++) {
            scanf("%d%d%d", &op, &x, &y);
            if (!op) {
                scanf("%d", &z);
                Min(1, x, y, z);
            } else {
                printf("%lld\n", query(1, x, y, op));
            }
        }
    }
    return 0;
}