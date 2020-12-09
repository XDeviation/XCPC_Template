#include <bits/stdc++.h>
#define T 100
#define eps 1e-8
#define delta 0.98
#define INF 0x7fffffff
using namespace std;
typedef long long ll;
const int maxn = 1e2 + 5;

struct point {
    double x, y, z;
} p[maxn];

double dis(point A, point B) {
    return sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y) +
                (A.z - B.z) * (A.z - B.z));
}

double GetSum(int n, point t) {
    double ans = 0;
    for (int i = 0; i < n; i++)
        ans += dis(t, p[i]);
    return ans;
}

double Search(int n) {
    point s = p[0];
    double t = T;
    double ans = INF;
    while (t > eps) {
        double maxx = 0;
        int k = 0;
        for (int i = 0; i < n; i++) {
            double temp = dis(p[i], s);
            if (temp > maxx) {
                maxx = temp;
                k = i;
            }
        }
        ans = min(ans, maxx);
        s.x += (p[k].x - s.x) / maxx * t;
        s.y += (p[k].y - s.y) / maxx * t;
        s.z += (p[k].z - s.z) / maxx * t;
        t *= delta;
    }
    return ans;
}

// 朴素的板子题 给你一堆点，求最小球覆盖

int main() {
    int n;
    while (scanf("%d", &n) && n) {
        for (int i = 0; i < n; i++)
            scanf("%lf%lf%lf", &p[i].x, &p[i].y, &p[i].z);
        double ans = Search(n);
        printf("%.5f\n", ans);
    }
    return 0;
}