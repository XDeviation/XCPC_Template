/*

   在这份代码实现中，rem表示的是len，0号节点的实际编号为1

   上文中提到的未插入后缀数rem没有被维护，而是被len和now表示

   len[u]和start[u]表示u的父边在s中的起点和长度，
   即u的父边代表s[start[u] ... start[u]+len[u]-1]

 */
#include <bits/stdc++.h>
using namespace std;
const int maxn = 1e5 + 7;
const int inf = -1e9;
struct suffixTree {
    int link[maxn << 1], len[maxn << 1], start[maxn << 1], s[maxn << 1], n,
        tail, now, rem, ch[maxn << 1][27];
    suffixTree()
        : tail(1)
        , n(0)
        , rem(0)
        , now(1) {
        len[0] = inf;
    }
    int newnode(int st, int le) {
        link[++tail] = 1;
        start[tail] = st;
        len[tail] = le;
        return tail;
    }
    void extend(int x) {
        s[++n] = x;
        rem++;
        for (int last = 1; rem;) {
            while (rem > len[ch[now][s[n - rem + 1]]])
                rem -= len[now = ch[now][s[n - rem + 1]]];
            int &v = ch[now][s[n - rem + 1]];
            int c = s[start[v] + rem - 1];
            if (!v || x == c) {
                link[last] = now;
                last = now;
                if (!v)
                    v = newnode(n, inf);
                else
                    break;
            } else {
                int u = newnode(start[v], rem - 1);
                ch[u][c] = v;
                ch[u][x] = newnode(n, inf);
                start[v] += rem - 1;
                len[v] -= rem - 1;
                link[last] = v = u;
                last = u;
            }
            if (now == 1)
                rem--;
            else
                now = link[now];
        }
    }
};
int main() {
}