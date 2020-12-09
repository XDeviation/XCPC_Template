#include <bits/stdc++.h>
using namespace std;
#define ll long long
const double eps = 1e-6;
const int maxn = 1e5 + 7;

// tips: this tmp maybe need each vector has intersection

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
} ans[maxn], tmpp[maxn];

double cross(point a, point b)
//求向量外积
{
    return a.x * b.y - a.y * b.x;
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
} li[maxn], tmpl[maxn];

point intersection(line a, line b)
// 交点
{
    return a.s + (a.t - a.s) * (cross(a.s - b.s, b.t - b.s) /
                                cross(b.t - b.s, a.t - a.s));
}

bool parallel(line a, line b)
//平行
{
    return dcmp(cross(a.t - a.s, b.t - b.s), 0);
}

bool isrig(line a, point b)
//判断点是否在向量右边
{
    return cross(a.t - a.s, b - a.s) + eps < 0;
}

bool linecmp(line a, line b)
// 极角序
{
    return dcmp(a.ang, b.ang) ? cross(a.t - a.s, b.t - a.s) + eps < 0
                              : a.ang < b.ang;
}

bool SI(line *li, int n, point *ret, int &m, line *ql, point *qp)
// 半平面交，li是向量集合，n是数量，ret是答案，ql和qp是临时的
{
    int l, r;
    sort(li + 1, li + 1 + n, linecmp);
    ql[l = r = 1] = li[1];
    for (int i = 2; i <= n; i++)
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

point a[maxn], b[maxn];

int main() {
    int t;
    cin >> t;
    int tot = 0;
    while (t--) {
        int n;
        cin >> n;
        for (int i = 0; i < n; i++)
            a[i].read();
        for (int i = 0; i < n; i++)
            li[++tot].getline(a[i], a[(i + 1) % n]);
    }
    int m = 0;
    if (SI(li, tot, ans, m, tmpl, tmpp))
        printf("%.3f", area(ans, m));
    else
        cout << "0.000";
    return 0;
}