#include <bits/stdc++.h>
using namespace std;
#define ll long long
const double eps = 1e-6;
const int maxn = 1e5 + 7;
#define sqr(x) ((x) * (x))
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

int pointcmp(const point &a, const point &b) {
    if (a.x != b.x)
        return (a.x < b.x);
    else
        return (a.y < b.y);
}

bool same_line(point x1, point x2, point x3)
//判断点是否在线段上，所求点为x1，若改为直线则在if中直接return true
{
    const double esp = 1e-8;
    if (abs((x3.y - x1.y) * (x2.x - x1.x) - (x2.y - x1.y) * (x3.x - x1.x)) <
        esp) {
        if (x1.x <= x2.x && x1.x >= x3.x && x1.y <= x2.y && x1.y >= x3.y)
            return true;
        if (x1.x >= x2.x && x1.x <= x3.x && x1.y >= x2.y && x1.y <= x3.y)
            return true;
        return false;
    }
    return false;
}

bool in_the_area(point p, point area[], int lens)
//判断点是否在多边形area内部，area内的点逆时针排序，lens为多边形点的数量，same_line用来区分边界点算内部还是外部
{
    int ans = 0;
    double x;
    for (int i = 1; i <= lens; i++) {
        point p1 = area[i];
        point p2 = area[i + 1];
        if (same_line(p, p1, p2)) return true;
        if (p1.y == p2.y) continue;
        if (p.y < min(p1.y, p2.y)) continue;
        if (p.y >= max(p1.y, p2.y)) continue;
        x = (p.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;
        if (x > p.x) ans++;
    }
    return (ans % 2 == 1);
}

ll multi(point p1, point p2, point p0)
// 叉乘
{
    ll x1 = p1.x - p0.x;
    ll y1 = p1.y - p0.y;

    ll x2 = p2.x - p0.x;
    ll y2 = p2.y - p0.y;

    return (x1 * y2 - x2 * y1);
}

int ConvexHull(int n, point p[], point ch[])
// 凸包，点p，输出ch，如果三点共线算凸包，那么把小于号改为小于等于号
{
    sort(p, p + n, pointcmp);
    int m = 0;
    for (int i = 0; i < n; i++) {
        while (m > 1 && multi(ch[m - 1], p[i], ch[m - 2]) < 0)
            m--;
        ch[m++] = p[i];
    }
    int k = m;
    for (int i = n - 2; i >= 0; i--) {
        while (m > k && multi(ch[m - 1], p[i], ch[m - 2]) < 0)
            m--;
        ch[m++] = p[i];
    }
    if (n > 1) m--;
    return m;
}

int cnt(int x, int y, int m) {
    x = x % m;
    y = y % m;
    if (x <= y)
        return y - x + 1;
    else
        return y + m - x + 1;
}

double Rotating_Caliper(int m, point ch[])
// 旋转卡壳，i，j是俩点
{
    ch[m + 1] = ch[1];
    int j = 2;
    double ans = 0;
    for (int i = 1; i <= m; i++) {
        while (fabs(multi(ch[i], ch[i + 1], ch[j])) <
               fabs(multi(ch[i], ch[i + 1], ch[j + 1]))) {
            j++;
            if (j > m) j = 1;
        }
        ans = max(ans, (ch[i] - ch[j]).len());
    }
    return ans;
}

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

double getd(point a, point b)
//求直径
{
    return sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));
}

point geto(point p1, point p2, point p3)
//求圆心
{
    point res;
    double a = p2.x - p1.x;
    double b = p2.y - p1.y;
    double c = p3.x - p2.x;
    double d = p3.y - p2.y;
    double e = sqr(p2.x) + sqr(p2.y) - sqr(p1.x) - sqr(p1.y);
    double f = sqr(p3.x) + sqr(p3.y) - sqr(p2.x) - sqr(p2.y);
    res.x = (f * b - e * d) / (c * b - a * d) / 2.0;
    res.y = (a * f - e * c) / (a * d - b * c) / 2.0;
    return res;
}

void mincir(point &o, double &r, point p[], int n)
// 最小圆覆盖
{
    o = p[1];
    r = 0;
    for (int i = 0; i < n; ++i) {
        if (getd(p[i], o) - r > eps) { //不在圆内
            o = p[i];
            r = 0;
            for (int j = 0; j < i; j++) {
                if (getd(p[j], o) - r > eps) { //不在圆内
                    o = (point){(p[i].x + p[j].x) / 2.0,
                                (p[i].y + p[j].y) / 2.0};
                    r = getd(p[i], p[j]) / 2.0;
                    for (int k = 0; k < j; ++k)
                        if (getd(p[k], o) - r > eps) {  //不在圆内
                            o = geto(p[i], p[j], p[k]); //外接圆
                            r = getd(p[i], o);
                        }
                }
            }
        }
    }
}

struct line
// 线
{
    point s, t;
    double a, b, c;
    double ang;
    void getline(point a, point b) {
        s = a, t = b;
        ang = (t - s).angle();
    }
    line() {
    }
    line(point x, point y) {
        this->s = x, this->t = y, this->ang = (x - y).angle();
        if (fabs(x.x - y.x) < eps) {
            a = 1, c = -x.x, b = 0;
        } else if (fabs(x.y - y.y) < eps) {
            a = 0, b = 1, c = -x.y;
        } else {
            a = x.y - y.y, b = y.x - x.x,
            c = (x.x - y.x) * y.y - (x.y - y.y) * y.x;
        }
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

int l, r;

bool SI(line *li, int n, point *ret, int &m, line *ql, point *qp)
// 半平面交，li是向量集合，n是数量，ret是答案，ql和qp是临时的，求的是左侧的半平面交
{
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

struct Point
// 三维
{
    double x, y, z;
    void read() {
        scanf("%lf%lf%lf", &x, &y, &z);
    }
    Point operator+(const Point &p) const {
        return {x + p.x, y + p.y, z + p.z};
    }
    Point operator-(const Point &p) const {
        return {x - p.x, y - p.y, z - p.z};
    }
    Point operator*(double p) const {
        return {x * p, y * p, z * p};
    }
    Point operator/(double p) const {
        return {x / p, y / p, z / p};
    }
    Point operator*(const Point &p) const {
        return {y * p.z - z * p.y, -x * p.z + z * p.x, x * p.y - y * p.x};
    }
    double operator^(const Point &p) const {
        return x * p.x + y * p.y + z * p.z;
    }
    double len() {
        return sqrt(x * x + y * y + z * z);
    }
    double dis(const Point &p) const {
        return ((*this) - p).len();
    }
};

point a[maxn], b[maxn];

int main() {
    return 0;
}