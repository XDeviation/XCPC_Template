#include <bits/stdc++.h>
#ifndef open_dbg_func
#define dbg(args...) (args)
#endif
using namespace std;
struct state {
    int len, link;
    map<char, int> nxt;
};
const int maxn = 1e6 + 7;
long long ans;
state st[maxn * 2];
int len[maxn * 2], siz[maxn * 2];
int sz, last;
vector<int> leni[maxn];
void sam_init() {
    for (int i = 0; i < sz; i++) {
        st[i].nxt.clear();
        siz[i] = 0;
    }
    st[0].len = 0;
    st[0].link = -1;
    sz = 1;
    last = 0;
}

void sam_extend(char c) {
    int cur = sz++;
    st[cur].len = st[last].len + 1;
    int p = last;
    while (p != -1 && !st[p].nxt.count(c)) {
        st[p].nxt[c] = cur;
        p = st[p].link;
    }
    siz[cur] = 1;
    if (p == -1) {
        st[cur].link = 0;
    } else {
        int q = st[p].nxt[c];
        if (st[p].len + 1 == st[q].len) {
            st[cur].link = q;
        } else {
            int clone = sz++;
            st[clone].len = st[p].len + 1;
            st[clone].nxt = st[q].nxt;
            st[clone].link = st[q].link;
            while (p != -1 && st[p].nxt[c] == q) {
                st[p].nxt[c] = clone;
                p = st[p].link;
            }
            st[q].link = st[cur].link = clone;
        }
    }
    last = cur;
}
int main() {
    string s;
    cin >> s;
    sam_init();
    for (int i = 0; i < s.length(); i++)
        sam_extend(s[i]);
    for (int i = 1; i < sz; i++) {
        leni[st[i].len].push_back(i);
    }
    ans = 0;
    for (int i = s.length(); i > 0; i--) {
        for (auto it : leni[i]) {
            siz[st[it].link] += siz[it];
            if (siz[it] > 1) ans = max(ans, 1LL * siz[it] * st[it].len);
        }
    }
    cout << ans;
}