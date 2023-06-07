#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

#define MatMax 20

//函数声明
double Det(const double arr[MatMax][MatMax], int n);
double Cof(const double arr[MatMax][MatMax], int i, int n);
void FindCof(double arr[MatMax][MatMax], double arr2[MatMax][MatMax], int i, int j, int n);
double **matrix_inver(double **arr, int n);
double **creat_metrix_n3(double src[3][3]);
void print_mat(double **arr, int row, int col);

//计算行列式
double Det(const double arr[MatMax][MatMax], int n)
{
    assert(n > 0);
    double sum = 0;
    int i = 0;
    if (n == 1) // 1阶行列式直接得出结果
    {
        sum = arr[0][0];
    }
    else if (n == 2)
    {
        sum = arr[0][0] * arr[1][1] - arr[0][1] * arr[1][0]; //杀戮法求解
    }
    else if (n == 3)
    {
        sum = arr[0][0] * arr[1][1] * arr[2][2] + arr[0][1] * arr[1][2] * arr[2][0] + arr[1][0] * arr[2][1] * arr[0][2] - arr[0][2] * arr[1][1] * arr[2][0] - arr[0][1] * arr[1][0] * arr[2][2] - arr[1][2] * arr[2][1] * arr[0][0]; //划线法求解
    }
    else
    {
        for (i = 0; i < n; i++) //按第一行展开
        {
            if (arr[0][i] != 0) //展开项不为0才计算
            {
                sum += ((int)pow(-1, i + 0)) * arr[0][i] * (Cof(arr, i, n)); // 2阶以上继续递归
            }
            else
                sum += 0; //展开项为0
        }
    }
    return sum;
}
//找到余子式
double Cof(const double arr[MatMax][MatMax], int i, int n)
{
    assert(n > 0);
    int k = 0;
    int j = 0;
    double arr2[MatMax][MatMax] = {0};
    for (k = 0; k < n - 1; k++) //去除0行i列，剩下的组成新的矩阵
    {
        for (j = 0; j < n - 1; j++)
        {
            if (j < i)
            {
                arr2[k][j] = arr[k + 1][j];
            }
            else
            {
                arr2[k][j] = arr[k + 1][j + 1];
            }
        }
    }
    return Det(arr2, n - 1);
}
//找到去掉i行j列的余子式
void FindCof(double arr[MatMax][MatMax], double arr2[MatMax][MatMax], int i, int j, int n)
{
    int m = 0;
    int k = 0;
    for (m = 0; m < n - 1; m++)
    {
        for (k = 0; k < n - 1; k++)
        {
            if (k < j)
            {
                if (m < i)
                {
                    arr2[m][k] = arr[m][k];
                }
                else
                {
                    arr2[m][k] = arr[m + 1][k];
                }
            }
            else
            {
                if (m < i)
                {
                    arr2[m][k] = arr[m][k + 1];
                }
                else
                {
                    arr2[m][k] = arr[m + 1][k + 1];
                }
            }
        }
    }
}

//计算逆的主函数
double **matrix_inver(double **arr, int n)
{
    int i, j;
    double **res = NULL;
    res = (double **)malloc(sizeof(double *) * n);
    if (res == NULL)
        exit(-1);
    for (i = 0; i < n; i++)
    {
        res[i] = (double *)malloc(sizeof(double) * n);
        memset(res[i], 0, sizeof(double) * n);
    }
    double tmp[MatMax][MatMax] = {0};
    //保护arr，将arr指向内存的数据拷贝到tmp二维数组中
    for (i = 0; i < n; i++)
    {
        memcpy(tmp[i], arr[i], sizeof(double) * n);
    }
    double a = 1.0 / (Det(tmp, n));
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            double tmp2[MatMax][MatMax] = {0};
            FindCof(tmp, tmp2, j, i, n); //求转置后的伴随
            double b = pow(-1, i + j) * Det(tmp2, n - 1);
            res[i][j] = a * b;
        }
    }
    return res;
}

//创建n维矩阵空间，并初始化
double **test1(int n)
{
    double **arr = (double **)malloc(sizeof(double *) * n);
    int i, j;
    if (arr != NULL)
    {
        for (i = 0; i < n; i++)
        {
            arr[i] = (double *)malloc(sizeof(double) * n);
        }
        //为矩阵赋值
        if (*arr != NULL)
        {
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    arr[i][j] = pow(i, j);
                }
            }
        }
    }
    return arr;
}

//创建n维矩阵空间，并初始化
double **creat_metrix_n3(double src[3][3])
{
    double **arr = (double **)malloc(sizeof(double *) * 3);
    int i, j;
    if (arr != NULL)
    {
        for (i = 0; i < 3; i++)
        {
            arr[i] = (double *)malloc(sizeof(double) * 3);
        }
        //为矩阵赋值
        if (*arr != NULL)
        {
            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    arr[i][j] = src[i][j];
                }
            }
        }
    }
    return arr;
}

//打印矩阵
void print_mat(double **arr, int row, int col)
{
    putchar('\n');
    int i, j;
    // row = (int)_msize(arr) / (int)sizeof(double*);//判断行数
    // col = (int)_msize(*arr) / (int)sizeof(double);//判断列数
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            printf("%10.5lf ", arr[i][j]);
        }
        putchar('\n');
    }
    putchar('\n');
}

void print_arr(double *arr, int len)
{
    putchar('\n');
    int i, j;
    // row = (int)_msize(arr) / (int)sizeof(double*);//判断行数
    // col = (int)_msize(*arr) / (int)sizeof(double);//判断列数
    for (i = 0; i < len; i++)
    {
        printf("%10.5lf ", arr[i]);
    }
    putchar('\n');
}

// int main()
// {
// 	int n = 5;
// 	double** arr = test1(n);
// 	printf("原矩阵:>\n");
// 	print(arr);
// 	double** res = matrix_inver(arr, n);
// 	printf("逆矩阵:>\n");
// 	print(res);
//     system("pause");

// 	return 0;
// }