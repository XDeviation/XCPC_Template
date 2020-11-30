#include <bits/stdc++.h>
using namespace std;
const int maxn = 505;
const double eps = 1e-8;
double A[maxn << 1][maxn],
    x[maxn]; // A矩阵中每一行1~n存系数，n+1为答案，m个方程m行，x是最终的答案
//注意空间要多开几个，还要考虑n，m不同的情况

int Guass(int n, int m) //有n个未知数，m个方程
{
    int i = 1, j = 1, k, r, c;
    while (i <= m && j <= n) //正在处理第i个方程，解第j个未知数
    {
        r = i; //找到绝对值最大的系数，防止除数为0的情况，使得其他方程组系数不会变得太大
        for (k = i + 1; k <= m; k++)
            if (fabs(A[k][j]) > fabs(A[r][j])) r = k;
        if (fabs(A[r][j]) >=
            eps) //出现为0的情况，说明此项已经被消掉了，直接用进行下一个未知数，而方程不变,不过这个时候，一般来说跳过的这个元素就没有固定解啦
        {
            for (c = 1; c <= n + 1; c++)
                swap(A[i][c], A[r][c]); //交换
            for (k = i + 1; k <= m; k++)
                if (fabs(A[k][j]) >= eps) {
                    double f = A[k][j] / A[i][j];
                    for (c = j; c <= n + 1; c++) //当前方程j前面的系数都是0
                        A[k][c] -= f * A[i][c];
                }
            i++; //获取下一个方程
        }
        j++; //去消下一个未知数
    }
    //必须先判无解再判断多解
    for (k = i; k <= m; k++)
        if (fabs(A[k][n + 1]) >= eps)
            return 0; //若有一行系数为0但是不为答案，则无解
    if (i <= n)
        return 2; //如果被你处理出来的方程没有n个，就会出现多解。(i=n表示解决了n-1个方程)
    for (int i = n; i >= 1; i--) {
        for (j = i + 1; j <= n; j++)
            A[i][n + 1] -= A[i][j] * x[j];
        x[i] = A[i][n + 1] / A[i][i];
    }
    //最终统计出来的答案x[i]肯定是对应的第i个元素的解哦,换的只是方程的顺序
    return 1; //拥有唯一解
}
// 模线性方程组 除法换成逆元

int B[maxn << 1][maxn], cnt[maxn], ans[maxn][maxn];
void xorGauss(int n, int m) {
    int i = 1, j = 1, k, r, c;
    while (i <= m && j <= n) // i是正在考虑的方程，j是待求解的系数
    {
        for (r = i; r <= m; r++)
            if (B[r][j]) break; //找到有值的数简单些
        if (r <= m)             //新判断
        {
            for (c = 1; c <= n + 1; c++)
                swap(B[r][c], B[i][c]);
            for (k = i + 1; k <= m; k++)
                if (B[k][j])
                    for (c = j; c <= n + 1; c++)
                        B[k][c] ^= B[i][c];
            i++;
        }
        j++;
    }
    i--, j--;
    //上面时解方程部分，下面是搜索自由元部分
    /*
    不搜索自由元可以这样
    for(k=i;k<=m;k++)if(B[k][n+1])return -1;//判断无解
    return n-i+1;//返回自由元个数
    */
    for (int k = 1; k <= i; k++)
        for (int j = n; j >= 1; j--)
            if (B[k][j]) cnt[k] = j;
    while (
        i >= 1 ||
        j >=
            1) //即使没有方程了，也要把自由元算完，但是不会出现还有方程却没有自由元的情况
    {
        if (cnt[i] != j) {
            ans[1][j] = 1;
            j--;
        } else {
            for (int k = j + 1; k <= n; k++)
                if (B[i][k]) B[i][n + 1] ^= ans[1][k];
            ans[1][j] = B[i][n + 1];
            i--, j--;
        }
    }
}
// 异或方程组

int main() {
    int n, m;
    cin >> n;
    m = n;
    for (int i = 1; i <= m; i++)
        for (int j = 1; j <= n + 1; j++)
            cin >> A[i][j];
    int ans = Guass(n, m);
    if (ans != 1)
        cout << "No Solution\n";
    else
        for (int i = 1; i <= n; i++)
            cout << fixed << setprecision(2) << x[i] << '\n';
}