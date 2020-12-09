#include <bits/stdc++.h>
#define For(i, a, b, c) for (int(i) = (a); (i) <= (b); (i) += (c))
#define Bor(i, a, b, c) for (int(i) = (b); (i) >= (a); (i) -= (c))
#define Min(a, b) ((a) > (b) ? (b) : (a))
#define Max(a, b) ((a) > (b) ? (a) : (b))
using namespace std;
const int N = 200050;
int p[N], ll[N], ans, rr[N], mx, id, cnt;
char s[N], t[N];
// 最长双回文串长度 双回文串就是两个回文串拼起来
int main() {
    scanf("%s", t);
    int len = strlen(t);
    s[++cnt] = '$', s[++cnt] = '#';
    For(i, 0, len - 1, 1) s[++cnt] = t[i], s[++cnt] = '#';
    s[++cnt] = '\0';
    For(i, 1, cnt, 1) {
        if (i < mx)
            p[i] = Min(p[id * 2 - i], mx - i);
        else
            p[i] = 1;
        while (s[i - p[i]] == s[i + p[i]])
            p[i]++;
        if (mx < i + p[i]) id = i, mx = i + p[i];
        ll[i + p[i] - 1] = Max(ll[i + p[i] - 1], p[i] - 1);
        rr[i - p[i] + 1] = Max(rr[i - p[i] + 1], p[i] - 1);
    }
    For(i, 2, cnt, 2) rr[i] = Max(rr[i], rr[i - 2] - 2);
    Bor(i, 2, cnt, 2) ll[i] = Max(ll[i], ll[i + 2] - 2);
    For(i, 2, cnt, 2) if (rr[i] && ll[i]) ans = Max(ans, ll[i] + rr[i]);
    cout << ans;
    return 0;
}