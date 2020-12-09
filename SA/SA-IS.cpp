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

const int N = 1e5 + 20;
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

// 读入一个长度为 n 的由小写英文字母组成的字符串，
// 请把这个字符串的所有非空后缀按字典序从小到大排序，
// 然后按顺序输出后缀的第一个字符在原串中的位置。位置编号为 1 到 n。
// 除此之外为了进一步证明你确实有给后缀排序的超能力，请另外输出
// n−1 个整数分别表示排序后相邻后缀的最长公共前缀的长度。

int main() {
    int n, i, j, h;
    n = fread(s, 1, N, stdin);
    while (s[n - 1] - 97u > 25)
        --n;
    for (i = 0; i < n; ++i)
        ord[i] = s[i] - 96;
    suffix_array(n, 26, ord, sa);
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
    for (i = 0; i < n; ++i) {
        Write::print(sa[i] + 1);
        Write::putchar(" \n"[i == n - 1]);
    }
    for (i = 1; i < n; ++i) {
        Write::print(height[i]);
        Write::putchar(" \n"[i == n - 1]);
    }
    Write::flush();
    return 0;
}