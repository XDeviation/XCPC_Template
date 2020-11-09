#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
const int BIT = 808, E = 62;
int n1, n2, m, h, i, j, p, ans;
int s[50000], t[50000];
struct Num {
    ll x[BIT];
    Num() {
        for (int i = 0; i < BIT; i++)
            x[i] = 0;
    }
    void set(int p) {
        x[p / E] |= 1LL << (p % E);
    }
    Num operator&(Num b) {
        Num c;
        for (int i = 0; i <= m; i++)
            c.x[i] = x[i] & b.x[i];
        return c;
    }
    Num operator|(Num b) {
        Num c;
        for (int i = 0; i <= m; i++)
            c.x[i] = x[i] | b.x[i];
        return c;
    }
    Num operator^(Num b) {
        Num c;
        for (int i = 0; i <= m; i++)
            c.x[i] = x[i] ^ b.x[i];
        return c;
    }
    Num operator-(Num b) {
        Num c;
        for (int i = 0; i <= m; i++)
            c.x[i] = x[i] - b.x[i];
        for (int i = 0; i < m; i++)
            if (c.x[i] < 0) c.x[i] += (1LL << E), c.x[i + 1];
        return c;
    }
    void shl() {
        ll o = 1, p;
        for (int i = 0; i <= m; o = p, i++) {
            p = x[i] & (1LL << h), (x[i] <<= 1) &= ~(1LL << (h + 1));
            if (o) x[i] |= 1;
        }
    }
} ap[100007], x, row[2];
int my_hash(int x) {
    if (x == 'A') return 0;
    if (x == 'C') return 1;
    if (x == 'G') return 2;
    return 3;
}
int main() {
    // scanf("%d%d%s%s", &n1, &n2, s, t);
    int n;
    cin >> n;
    n1 = n, n2 = n;
    for (int it = 1; it <= n; it++)
        cin >> s[it];
    for (int it = 1; it <= n; it++)
        cin >> t[it];
    i = 0;
    for (m = (n1 - 1) / E, h = (m ? E : n1) - 1; i < n1; i++)
        ap[s[i]].set(i);
    for (i = 0; i < n2; i++) {
        p ^= 1;
        x = ap[t[i]] | row[p ^ 1];
        row[p ^ 1].shl();
        row[p] = x & ((x - row[p ^ 1]) ^ x);
    }
    for (i = m; ~i; i--)
        for (j = h; ~j; j--)
            if (row[p].x[i] & (1LL << j)) ans++;
    printf("%d", ans);
}