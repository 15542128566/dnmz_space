#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include "common.c"

// typedef struct mat
// {
//     int rows;
//     int cols;
//     unsigned char **data;
// } Mat;

// 定义矩阵的行数和列数
#define ROWS 3
#define COLS 3

// Mat create_mat(int rows, int cols)
// {
//     Mat mat;
//     mat.data = NULL;
//     mat.rows = rows;
//     mat.cols = cols;
//     mat.data = (unsigned char **)malloc(rows * sizeof(unsigned char *));
//     for (int i = 0; i < rows; ++i)
//     {
//         mat.data[i] = (unsigned char *)malloc(cols * sizeof(unsigned char));
//         memset(mat.data[i], 0, cols * sizeof(unsigned char));
//     }
//     return mat;
// }

int main()
{
    // 构建 Mat 的数组
    int num_mats = 2;
    Mat mats[num_mats];
    for (int i = 0; i < num_mats; i++)
    {
        mats[i] = create_mat(ROWS, COLS);
        for (int j = 0; j < ROWS; j++)
        {
            for (int k = 0; k < COLS; k++)
            {
                mats[i].data[j][k] = 1;
            }
        }
    }

    // 输出 Mat 数组的值
    for (int i = 0; i < num_mats; i++)
    {
        printf("Mat %d:\n", i);
        for (int j = 0; j < mats[i].rows; j++)
        {
            for (int k = 0; k < mats[i].cols; k++)
            {
                printf("%d ", mats[i].data[j][k]);
            }
            printf("\n");
        }
    }

    // 释放 Mat 数组的空间
    for (int i = 0; i < num_mats; i++)
    {
        for (int j = 0; j < mats[i].rows; j++)
        {
            free(mats[i].data[j]);
        }
        free(mats[i].data);
    }

    system("pause");

    return 0;
}
