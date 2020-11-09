/*
 * Codeforces Contest
 * Au: Lost_Deviation
 * Time: 2020-11-09 15:15:05
 */
#include "dbg_func"
#include <bits/extc++.h>
#include <bits/stdc++.h>
/*
若编译器没有 extc++.h 使用如下三个库
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/priority_queue.hpp>
#include <ext/pb_ds/tree_policy.hpp>
*/
#define param int, less<int>, pairing_heap_tag
using namespace std;
// using namespace __gnu_pbds;
#ifndef open_dbg_func
#define dbg(args...) (args)
#endif
using namespace std;
#define ll long long
#define ull unsigned long long
const int maxn = 1e6 + 7;
const double eps = 1e-6;
const double pi = acos(-1);
const ll p = 998244353;
const int maxp = 300010;
static auto fast_io = []() {
    std::ios::sync_with_stdio(false); // turn off sync
    cin.tie(nullptr);                 // untie in/out streams
    return 0;
}();
vector<string> split(string str, const string &delimiters) {
    regex del(delimiters);
    vector<string> v(sregex_token_iterator(str.begin(), str.end(), del, -1),
                     sregex_token_iterator());
    return v;
}
bool dcmp(double a, double b) {
    return fabs(a - b) < eps;
}
class point
// 点或向量的类
{
  public:
    double x, y;
    void read() {
        cin >> x >> y;
    }
    void print() {
        cout << x << y;
    }

    point operator+(const point &p) const {
        return {x + p.x, y + p.y};
    }
    point operator-(const point &p) const {
        return {x - p.x, y - p.y};
    }
    point operator*(double p) const {
        return {x * p, y * p};
    }
    point operator/(double p) const {
        return {x / p, y / p};
    }
    double operator*(const point &p) const {
        return x * p.x + y * p.y;
    }
    double len() {
        return sqrt(x * x + y * y);
    }
    double dis(const point &p) const {
        return ((*this) - p).len();
    }
    double angle() {
        return atan2(y, x);
    }
} a[maxn], tmpp[maxn], ans[maxn];
struct line
// 线
{
    point s, t;
    double ang;
    void getline(point a, point b) {
        s = a, t = b;
        ang = (t - s).angle();
    }
    line() {
    }
    line(point a, point b) {
        this->s = a, this->t = b, this->ang = (a - b).angle();
    }
    double len() {
        return (s - t).len();
    }
} li[maxn], tmpl[maxn];
double drgcos(line a, line b) {
    return ((b.t - b.s) * (a.t - a.s) / a.len() * b.len());
}
double cross(point a, point b)
//求向量外积
{
    return a.x * b.y - a.y * b.x;
}
bool parallel(line a, line b)
//平行
{
    return dcmp(cross(a.t - a.s, b.t - b.s), 0);
}
bool linecmp(line a, line b)
// 极角序
{
    return dcmp(a.ang, b.ang) ? cross(a.t - a.s, b.t - a.s) + eps < 0
                              : a.ang < b.ang;
}
bool isrig(line a, point b)
//判断点是否在向量右边
{
    return cross(a.t - a.s, b - a.s) + eps < 0;
}
point intersection(line a, line b)
// 交点
{
    return a.s + (a.t - a.s) * (cross(a.s - b.s, b.t - b.s) /
                                cross(b.t - b.s, a.t - a.s));
}
bool SI(line *li, int n, point *ret, int &m, line *ql, point *qp)
// 半平面交，li是向量集合，n是数量，ret是答案，ql和qp是临时的
{
    int l, r;
    sort(li + 1, li + 1 + n, linecmp);
    ql[l = r = 1] = li[1];
    for (int i = 2; i <= n; i++) {
        dbg(i);
        if (!dcmp(li[i].ang, li[i - 1].ang)) {
            while (l < r && isrig(li[i], qp[r - 1]))
                --r;
            while (l < r && isrig(li[i], qp[l]))
                ++l;
            if (l <= r) qp[r] = intersection(ql[r], li[i]);
            ql[++r] = li[i];
            if (l < r &&
                (parallel(ql[l], ql[l + 1]) || parallel(ql[r], ql[r - 1])))
                return false;
        }
    }
    // dbg(n);
    while (l < r && isrig(ql[l], qp[r - 1]))
        --r;
    while (l < r && isrig(ql[r], qp[l]))
        ++l;
    if (r - l <= 1) return false;
    qp[r] = intersection(ql[l], ql[r]);
    m = 0;
    for (int i = l; i <= r; i++)
        ret[++m] = qp[i];
    return true;
}
double solvef(double a, double b, double c) {
    return max((-b + sqrt(b * b - 4 * a * c)) / (2 * a),
               (-b - sqrt(b * b - 4 * a * c)) / (2 * a));
}
double area(point s[], int n)
//求凸多边形面积
{
    double ret = 0;
    s[n + 1] = s[1];
    for (int i = 1; i <= n; i++)
        ret += cross(s[i], s[i + 1]);
    return fabs(ret / 2);
}
signed main(void) {
    int n;
    cin >> n;
    for (int i = 0; i < n; i++) {
        a[i].read();
    }
    line str(a[0], a[1]);
    int tot = 0;
    for (int i = 1; i < n; i++) {
        line tmp(a[i], a[(i + 1) % n]);
        if (parallel(str, tmp)) {
            point inter =
                (str.t + tmp.s) * (tmp.len()) / (tmp.len() + str.len());
            point ender =
                (str.s + tmp.t) * (tmp.len()) / (tmp.len() + str.len());
            li[++tot] = line(inter, ender);
        } else {
            point inter = intersection(str, tmp);
            double x = cos(pi - acos(drgcos(str, tmp))),
                   k = tmp.len() / str.len();

            double y = (1 - x * x) / (2 * x + k * k + 1);
            // dbg(x, k, inter.x, inter.y, y);

            point add;
            add.x = 1, add.y = 1 / tan(tmp.ang + asin(sqrt(y)));
            li[++tot] = line(inter, inter - add);
        }
        dbg(tot, li[tot].s.x, li[tot].s.y, li[tot].t.x, li[tot].t.y,
            li[tot].ang);
    }
    int m;
    for (int i = 1; i <= tot; i++) {
    }
    SI(li, tot, ans, m, tmpl, tmpp);
    dbg(m);
    cout << fixed << setprecision(10) << area(ans, m) / area(a, n);
}
