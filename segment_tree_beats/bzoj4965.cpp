#include <bits/stdc++.h>
using namespace std;
#define N 500010
#define LL long long
#define INF 1 << 30
int n, m, a[N];
#define lc (p << 1)
#define rc (p << 1 | 1)
struct node {
    int l, r, len, mx, mn, sx, sn, cx, cn;
    LL laz, sum;
} t[N << 2];
void up(int p) {
    t[p].cx = t[p].cn = 0;
    t[p].sum = t[lc].sum + t[rc].sum;
    t[p].mx = max(t[lc].mx, t[rc].mx);
    t[p].sx = max(t[lc].sx, t[rc].sx);
    if (t[lc].mx ^ t[rc].mx) t[p].sx = max(t[p].sx, min(t[lc].mx, t[rc].mx));
    if (t[p].mx == t[lc].mx) t[p].cx += t[lc].cx;
    if (t[p].mx == t[rc].mx) t[p].cx += t[rc].cx;
    t[p].mn = min(t[lc].mn, t[rc].mn);
    t[p].sn = min(t[lc].sn, t[rc].sn);
    if (t[lc].mn ^ t[rc].mn) t[p].sn = min(t[p].sn, max(t[lc].mn, t[rc].mn));
    if (t[p].mn == t[lc].mn) t[p].cn += t[lc].cn;
    if (t[p].mn == t[rc].mn) t[p].cn += t[rc].cn;
}
void build(int p, int l, int r) {
    t[p].l = l;
    t[p].r = r;
    t[p].len = r - l + 1;
    if (l == r) {
        t[p].mx = t[p].mn = t[p].sum = a[l];
        t[p].cx = t[p].cn = 1;
        t[p].sx = -INF;
        t[p].sn = INF;
        return;
    }
    int mid = (l + r) >> 1;
    build(lc, l, mid);
    build(rc, mid + 1, r);
    up(p);
}
void now(int p, int v) {
    t[p].sum += 1ll * t[p].len * v;
    t[p].laz += v;
    t[p].mx += v;
    t[p].mn += v;
    if (t[p].sx != -INF) t[p].sx += v;
    if (t[p].sn != INF) t[p].sn += v;
}
void MAX(int p, int v) {
    t[p].sum += 1ll * (v - t[p].mn) * t[p].cn;
    t[p].mn = v;
    t[p].mx = max(v, t[p].mx);
    if (t[p].mx == t[p].mn) {
        t[p].sum = 1ll * t[p].len * v;
        t[p].cn = t[p].cx = t[p].len;
        t[p].sx = -INF;
        t[p].sn = INF;
    } else
        t[p].sx = max(v, t[p].sx);
}
void MIN(int p, int v) {
    t[p].sum += 1ll * (v - t[p].mx) * t[p].cx;
    t[p].mx = v;
    t[p].mn = min(v, t[p].mn);
    if (t[p].mx == t[p].mn) {
        t[p].sum = 1ll * t[p].len * v;
        t[p].cn = t[p].cx = t[p].len;
        t[p].sx = -INF;
        t[p].sn = INF;
    } else
        t[p].sn = min(v, t[p].sn);
}
void down(int p) {
    if (t[p].laz) {
        now(lc, t[p].laz);
        now(rc, t[p].laz);
        t[p].laz = 0;
    }
    if (t[lc].mn < t[p].mn && t[p].mn < t[lc].sn) MAX(lc, t[p].mn);
    if (t[rc].mn < t[p].mn && t[p].mn < t[rc].sn) MAX(rc, t[p].mn);
    if (t[lc].sx < t[p].mx && t[p].mx < t[lc].mx) MIN(lc, t[p].mx);
    if (t[rc].sx < t[p].mx && t[p].mx < t[rc].mx) MIN(rc, t[p].mx);
}
void add(int p, int l, int r, int v) {
    if (!v) return;
    if (l <= t[p].l && t[p].r <= r) {
        now(p, v);
        return;
    }
    int mid = (t[p].l + t[p].r) >> 1;
    down(p);
    if (l <= mid) add(lc, l, r, v);
    if (r > mid) add(rc, l, r, v);
    up(p);
}
void Max(int p, int l, int r, int v) {
    if (t[p].mn >= v) return;
    if (l <= t[p].l && t[p].r <= r && v < t[p].sn) {
        MAX(p, v);
        return;
    }
    int mid = (t[p].l + t[p].r) >> 1;
    down(p);
    if (l <= mid) Max(lc, l, r, v);
    if (r > mid) Max(rc, l, r, v);
    up(p);
}
void Min(int p, int l, int r, int v) {
    if (t[p].mx <= v) return;
    if (l <= t[p].l && t[p].r <= r && v > t[p].sx) {
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
        if (op == 4) return t[p].sum;
        if (op == 5) return t[p].mx;
        if (op == 6) return t[p].mn;
    }
    int mid = (t[p].l + t[p].r) >> 1;
    LL ans;
    down(p);
    if (op == 4) {
        ans = 0;
        if (l <= mid) ans += query(lc, l, r, op);
        if (r > mid) ans += query(rc, l, r, op);
    }
    if (op == 5) {
        ans = -INF;
        if (l <= mid) ans = max(ans, query(lc, l, r, op));
        if (r > mid) ans = max(ans, query(rc, l, r, op));
    }
    if (op == 6) {
        ans = INF;
        if (l <= mid) ans = min(ans, query(lc, l, r, op));
        if (r > mid) ans = min(ans, query(rc, l, r, op));
    }
    up(p);
    return ans;
}
/*
    1.给一个区间 [L, R] 加上一个数x
　　2.把一个区间 [L, R] 里小于 x 的数变成x
　　3.把一个区间 [L, R] 里大于 x 的数变成x
　　4.求区间 [L, R] 的和
　　5.求区间 [L, R] 的最大值
　　6.求区间 [L, R] 的最小值
*/
int main() {
    int op, x, y, z;
    scanf("%d", &n);
    for (int i = 1; i <= n; i++)
        scanf("%d", &a[i]);
    scanf("%d", &m);
    build(1, 1, n);
    for (int i = 1; i <= m; i++) {
        scanf("%d%d%d", &op, &x, &y);
        if (op <= 3) {
            scanf("%d", &z);
            if (op == 1) add(1, x, y, z);
            if (op == 2) Max(1, x, y, z);
            if (op == 3) Min(1, x, y, z);
        } else {
            printf("%lld\n", query(1, x, y, op));
        }
    }
    return 0;
}