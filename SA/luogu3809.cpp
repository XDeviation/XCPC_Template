#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define meow(args...) fprintf(stderr, args)
namespace Write {
const int N = 1 << 25;
char buf[N], s[20], *w = buf;
void flush() {
    fwrite(buf, 1, w - buf, stdout);
    w = buf;
}
inline void putchar(register char c) {
    *w++ = c;
}
void print(register int n) {
    register char *t = s;
    do {
        *t++ = n % 10 + 48;
    } while (n /= 10);
    while (t-- > s)
        putchar(*t);
}
} // namespace Write

const int N = 1e6 + 20;
int ord[2 * N], sa[2 * N], rank[N], height[N], l[N], r[N], pool4[4 * N],
    *mem4 = pool4;
char s[N], pool1[2 * N], *mem1 = pool1;
void induce(int n, int m, int *s, char *type, int *sa, int *cnt) {
    int i;
    memcpy(l + 1, cnt, m * sizeof(int));
    memcpy(r + 1, cnt + 1, m * sizeof(int));
    sa[l[s[n - 1]]++] = n - 1;
    for (i = 0; i < n; ++i) {
        int t = sa[i] - 1;
        if (t >= 0 && type[t]) sa[l[s[t]]++] = t;
    }
    for (i = n; i--;) {
        int t = sa[i] - 1;
        if (t >= 0 && !type[t]) sa[--r[s[t]]] = t;
    }
}
void suffix_array(int n, int m, int *s, int *sa) {
    int *cnt, *lms, *s1 = s + n, *sa1 = sa + n, n1 = 0, m1 = 0, i, t;
    char *type;
    type = mem1;
    mem1 += n + 1;
    cnt = mem4;
    mem4 += m + 1;
    type[n] = false;
    for (i = n; i--;) {
        type[i] = s[i] > s[i + 1] || (s[i] == s[i + 1] && type[i + 1]);
        ++cnt[s[i]];
    }
    for (i = 1; i <= m; ++i)
        r[i] = cnt[i] += cnt[i - 1];
    memset(rank, -1, n * sizeof(int));
    lms = mem4;
    for (i = 0; i < n; ++i) {
        if (!type[i] && (i == 0 || type[i - 1])) lms[rank[i] = n1++] = i;
    }
    lms[n1] = n;
    mem4 += n1 + 1;
    memset(sa, -1, n * sizeof(int));
    for (i = 0; i < n1; ++i)
        sa[--r[s[lms[i]]]] = lms[i];
    induce(n, m, s, type, sa, cnt);
    for (i = 0, t = -1; i < n; ++i) {
        int r = rank[sa[i]];
        if (r != -1) {
            int len = lms[r + 1] - sa[i] + 1;
            m1 += t == -1 || len != lms[rank[t] + 1] - t + 1 ||
                  memcmp(s + t, s + sa[i], len * sizeof(int));
            s1[r] = m1;
            t = sa[i];
        }
    }
    if (n1 == m1) {
        for (i = 0; i < n1; ++i)
            sa1[s1[i] - 1] = i;
    } else
        suffix_array(n1, m1, s1, sa1);
    memset(sa, -1, n * sizeof(int));
    memcpy(r + 1, cnt + 1, m * sizeof(int));
    for (i = n1; i--;) {
        t = lms[sa1[i]];
        sa[--r[s[t]]] = t;
    }
    induce(n, m, s, type, sa, cnt);
}
int main() {
    int n, i, j, h;
    n = fread(s, 1, N, stdin);
    while (s[n - 1] > 122u)
        --n;
    for (i = 0; i < n; ++i) {
        if (s[i] <= '9')
            ord[i] = s[i] - 47;
        else if (s[i] <= 'Z')
            ord[i] = s[i] - 65 + 11;
        else
            ord[i] = s[i] - 97 + 37;
    }
    suffix_array(n, 26 + 26 + 10, ord, sa);
    for (i = 0; i < n; ++i)
        rank[sa[i]] = i;
    for (i = 0, h = 0; i < n; ++i) {
        if (rank[i]) {
            j = sa[rank[i] - 1];
            while (s[i + h] == s[j + h])
                ++h;
            height[rank[i]] = h;
        } else
            h = 0;
        if (h) --h;
    }
    // Write::print(sa[4] + 1);
    for (i = 0; i < n; ++i) {
        Write::print(sa[i] + 1);
        Write::putchar(" \n"[i == n - 1]);
    }
    Write::flush();
    return 0;
}