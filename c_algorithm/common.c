#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// 宏定义
#define INPUT_FILE "test.bmp"
// #define INPUT_FILE "output0.bmp"
#define RGB_COUNT 3            // 每个像素包含色彩数
#define BMP_FILE_HEADER_LEN 14 // 位图文件BMP格式，文件头14字节
#define BMP_INFO_HEADER_LEN 40 // 位图文件BMP格式，信息头40字节

static const float COF_R2GRAY = 0.299, COF_G2GRAY = 0.587, COF_B2GRAY = 0.114; // rgb转化为灰度像素的系数

typedef struct mat
{
    int rows;
    int cols;
    unsigned char **data;
} Mat;

Mat create_mat(int rows, int cols)
{
    Mat mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (unsigned char **)malloc(rows * sizeof(unsigned char *));
    for (int i = 0; i < rows; ++i)
    {
        mat.data[i] = (unsigned char *)malloc(cols * sizeof(unsigned char));
        memset(mat.data[i], 0, cols * sizeof(unsigned char));
    }
    return mat;
}

Mat *create_mat_ptr(int rows, int cols)
{
    Mat *mat = (Mat *)malloc(sizeof(Mat));
    mat->rows = rows;
    mat->cols = cols;
    mat->data = (unsigned char **)malloc(rows * sizeof(unsigned char *));
    for (int i = 0; i < rows; ++i)
    {
        mat->data[i] = (unsigned char *)malloc(cols * sizeof(unsigned char));
        memset(mat->data[i], 0, cols * sizeof(unsigned char));
    }
    return mat;
}

void release_mat(Mat *mat)
{
    if (mat)
    {
        for (int i = 0; i < mat->rows; i++)
        {
            free(mat->data[i]);
            mat->data[i] = NULL;
        }
        free(mat->data);
        mat->data = NULL;
        free(mat);
        mat = NULL;
    }
}

typedef struct double_mat
{
    int rows;
    int cols;
    double **data;
} DoubleMat;

DoubleMat *create_double_mat_ptr(int rows, int cols)
{
    DoubleMat *double_mat = (DoubleMat *)malloc(sizeof(DoubleMat));
    double_mat->rows = rows;
    double_mat->cols = cols;
    double_mat->data = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; ++i)
    {
        double_mat->data[i] = (double *)malloc(cols * sizeof(double));
        memset(double_mat->data[i], 0, cols * sizeof(double));
    }
    return double_mat;
}
DoubleMat create_double_mat(int rows, int cols)
{
    DoubleMat double_mat;
    double_mat.rows = rows;
    double_mat.cols = cols;
    double_mat.data = (double **)malloc(rows * sizeof(double));
    for (int i = 0; i < rows; ++i)
    {
        double_mat.data[i] = (double *)malloc(cols * sizeof(double));
        memset(double_mat.data[i], 0, cols * sizeof(double));
    }
    return double_mat;
}

void release_double_mat_ptr(DoubleMat *double_mat)
{
    if (double_mat)
    {
        for (int i = 0; i < double_mat->rows; i++)
        {
            free((*double_mat).data[i]);
            (*double_mat).data[i] = NULL;
        }
        free((*double_mat).data);
        (*double_mat).data = NULL;
        free(double_mat);
        double_mat = NULL;
    }
}

Mat **convert_double_mat_to_mat_ptr2(DoubleMat **src, int rows, int cols)
{
    Mat **dest = (Mat **)malloc(rows * sizeof(Mat *));
    for (int i = 0; i < rows; ++i)
    {
        dest[i] = (Mat *)malloc(cols * sizeof(Mat));
        for (int j = 0; j < cols; ++j)
        {
            dest[i][j].rows = src[i][j].rows;
            dest[i][j].cols = src[i][j].cols;
            dest[i][j] = create_mat(src[i][j].rows, src[i][j].cols);
            for (int k = 0; k < src[i][j].rows; ++k)
            {
                for (int l = 0; l < src[i][j].cols; ++l)
                {
                    dest[i][j].data[k][l] = (unsigned char)src[i][j].data[k][l];
                }
            }
        }
    }
    return dest;
}

Mat *convert_double_mat_to_mat(const DoubleMat *src)
{

    Mat *dest = create_mat_ptr(src->rows, src->cols);
    for (int i = 0; i < src->rows; i++)
    {
        for (int j = 0; j < src->cols; j++)
        {
            dest->data[i][j] = (*src).data[i][j];
        }
    }
    return dest;
}

// Mat** new_mat_array(int lens) {
//     Mat** array = (Mat**) malloc(sizeof(Mat*) * lens);
//     for (int i = 0; i < rows; i++) {
//         array[i] = creat_mat(1, 1);
//     }
//     return array;
// }

// void delete_mat_array(Mat** array, int rows) {
//     for (int i = 0; i < rows; i++) {
//         delete_mat(array[i]);
//     }
//     free(array);
// }

// rgbToGray将图像处理为灰度图像
// 共col*row个像素，每个像素包含rgb3个字节，逐行逐列个像素处理
unsigned char *rgbToGray(unsigned char *src, int col, int row)
{
    unsigned char *gray = (unsigned char *)malloc(col * row * sizeof(unsigned char));
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            int idx = (i * col + j) * RGB_COUNT;
            int r = src[idx];
            int g = src[idx + 1];
            int b = src[idx + 2];
            gray[i * col + j] = COF_R2GRAY * r + COF_G2GRAY * g + COF_B2GRAY * b;
        }
    }

    return gray;
}

// 降采样函数
void downsample(unsigned char *img, int w, int h, int step, Mat *out)
{
    int new_w = w / step, new_h = h / step;
    // *out = create_mat(new_w, new_h);
    for (int i = 0; i < new_h; i++)
    {
        for (int j = 0; j < new_w; j++)
        {
            (*out).data[i][j] = img[i * w * step + j * step];
        }
    }

    // return out;
}

// 高斯核
void guassian_kernel(double sigma, int dim, DoubleMat *kernel)
{
    // Mat kernel = create_mat(dim, dim);
    (*kernel).rows = dim;
    (*kernel).cols = dim;
    double sum = 0.0;
    int r = dim / 2;
    for (int i = -r; i <= r; i++)
    {
        for (int j = -r; j <= r; j++)
        {
            double x = i * i + j * j;
            double value = exp(-x / (2 * sigma * sigma)) / (2 * M_PI * sigma * sigma);
            (*kernel).data[i + r][j + r] = value;
            // printf("kernel-->i=%d,j=%d,value=%f\n", i + r, j + r, (*kernel).data[i + r][j + r]);
            sum += value;
        }
    }
    // 归一化
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            // printf("before_kernel-->i=%d,j=%d,value=%f\n", i, j, (*kernel).data[i][j]);
            // printf("sum-->i=%d,j=%d,value=%f\n", i, j, sum);
            (*kernel).data[i][j] /= sum;
            // printf("normalized_kernel-->i=%d,j=%d,value=%f\n", i, j, (*kernel).data[i][j]);
        }
    }
}

// // 卷积计算
// // 卷积核中心点默认在该卷积核的中央，anchorx,anchory 取 0
// void convolve(Mat *src, DoubleMat *kernel, Mat *dest, int anchorx, int anchory)
// {
//     int r = kernel->rows / 2;
//     // *dest = create_mat(rows, cols);
//     for (int i = r; i < src->rows - r; i++)
//     {
//         for (int j = r; j < src->cols - r; j++)
//         {
//             int sum = 0;
//             for (int k = -r; k < r; k++)
//             {
//                 for (int l = -r; l < r; l++)
//                 {
//                     sum += kernel->data[k + r][l + r] * src->data[i + k - anchorx][j + l - anchory];
//                 }
//             }
//             dest->data[i][j] = sum;
//         }
//     }
// }

// 卷积计算
// 卷积核中心点默认在该卷积核的中央，anchorx,anchory 取 0
void convolve(Mat *src, DoubleMat *kernel, DoubleMat *dest, int anchorx, int anchory)
{
    int r = kernel->rows / 2;
    // *dest = create_mat(rows, cols);
    for (int i = 0; i < src->rows; i++)
    {
        for (int j = 0; j < src->cols; j++)
        {
            double sum = 0;
            for (int k = -r; k < r; k++)
            {
                for (int l = -r; l < r; l++)
                {
                    if (k + i >= 0 && k + i < src->rows && l + j >= 0 && l + j < src->cols)
                    {
                        sum += kernel->data[k + r][l + r] * src->data[i + k - anchorx][j + l - anchory];
                    }
                }
            }
            dest->data[i][j] = sum;
        }
    }
}

unsigned char *read_rgb_bmp(int *width, int *height)
{
    // 1.读取图像
    FILE *file = fopen(INPUT_FILE, "rb"); // rb-以二进制模式读取
    if (file == NULL)
    {
        printf("[main err] source image file null\n");
        return NULL;
    }
    // 获取图像大小
    int col = 0, row = 0;                           // 图像宽高
    fseek(file, BMP_FILE_HEADER_LEN + 4, SEEK_SET); // 跳过文件头和信息头的第一个字段
    fread(&col, sizeof(int), 1, file);
    fread(&row, sizeof(int), 1, file);
    // 存储图像数据
    unsigned char *src = NULL;
    src = (unsigned char *)malloc(col * row * RGB_COUNT * sizeof(unsigned char));
    fseek(file, (BMP_FILE_HEADER_LEN + BMP_INFO_HEADER_LEN), SEEK_SET); // 跳过文件头和信息头
    // todo: 应该根据颜色表判断的存储方式为RGB还是索引方式，这里简化处理为默认RGB存储方式
    fread(src, sizeof(char), col * row * RGB_COUNT, file);
    printf("width: %d, height: %d\n", col, row);
    *width = col;
    *height = row;
    fclose(file);

    return src;
}

// 四舍五入
int custom_round(double number)
{
    int rounded_num = (int)(round(number * 10 + 5) / 10);
    return rounded_num;
}

// 求最大值
double get_double_max(const double *src, int length)
{
    double res = src[0];
    for (int i = 1; i < length; i++)
    {
        if (src[i] > res)
        {
            res = src[i];
        }
    }

    return res;
}