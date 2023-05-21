#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define RGB_COUNT 3            // 每个像素包含色彩数
#define PIXEL_VALUE_MIN 0      // 像素最小范围
#define PIXEL_VALUE_MAX 255    // 像素最大范围
#define BMP_FILE_HEADER_LEN 14 // 位图文件BMP格式，文件头14字节
#define BMP_INFO_HEADER_LEN 40 // 位图文件BMP格式，信息头40字节
#define PI 3.14

static const float COF_R2GRAY = 0.299, COF_G2GRAY = 0.587, COF_B2GRAY = 0.114; // rgb转化为灰度像素的系数

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

// normalizeImage归一化处理，将图像中的每个像素缩放至[0~255]范围
unsigned char *normalizeImage(unsigned char *src, int col, int row)
{
    unsigned char *dst = (unsigned char *)malloc(col * row * sizeof(unsigned char));
    // 找到所有像素点中的最大最小值, 计算缩放比例scale
    int min = PIXEL_VALUE_MAX;
    int max = PIXEL_VALUE_MIN;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            int pv = src[i * col + j];
            min = (pv < min) ? pv : min;
            max = (pv > max) ? pv : max;
        }
    }
    if (max > min)
    {
        double scale = PIXEL_VALUE_MAX / (max - min);
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                dst[i * col + j] = (unsigned char)((src[i * col + j] - min) * scale);
            }
        }
    }
    else
    {
        printf("[normalizeImage warn] max value < min value");
        memcpy(dst, src, col * row); // 若max<min则放弃归一化，直接返回原始值
    }

    return dst;
}

// 均值滤波(边缘做截断处理) winsize为窗口边数
unsigned char *meanFilter(unsigned char *src, int col, int row, int winsize)
{
    if (col < winsize || row < winsize)
    {
        printf("averageFilter error: invalid params, col:%d,row:%d,winsize:%d\n", col, row, winsize);

        return NULL;
    }
    if (~(winsize & 1))
    {
        winsize += 1; //窗口边数调整为奇数
    }
    int winsize_2 = winsize / 2;
    unsigned char *dst = (unsigned char *)malloc(col * row * sizeof(unsigned char));
    memcpy(dst, src, col * row); // 边界值不会更改，copy原值
    double sum;
    int idx;
    for (int i = winsize_2; i < row - winsize_2; i++)
    {
        for (int j = winsize_2; j < col - winsize_2; j++)
        {
            sum = 0.0;
            for (int k = 0; k < winsize; k++)
            {
                for (int s = 0; s < winsize; s++)
                {
                    idx = (k + i - winsize_2) * col + j + s - winsize_2;
                    sum += src[idx];
                }
            }
            dst[i * col + j] = (unsigned char)(sum / (winsize * winsize));
        }
    }

    return dst;
}

// 边缘检测：基于sobel算子计算像素的梯度 kernel_gx = [-1 0 1 -2 0 2 -1 0 1] kernel_gy = [-1 -2 -1 0 0 0 1 2 1]
unsigned char *sobelDetectEgde(unsigned char *src, int col, int row)
{
    unsigned char *dst = (unsigned char *)malloc(col * row * sizeof(unsigned char));
    memset(dst, 0, col * row * sizeof(unsigned char));
    int gx, gy;
    for (int i = 1; i < row - 1; i++)
    {
        for (int j = 1; j < col - 1; j++)
        {
            // 计算水平方向梯度
            gx = -src[(i - 1) * col + j - 1] + src[(i - 1) * col + j + 1] - 2 * src[i * col + j - 1] +
                 2 * src[i * col + j + 1] - src[(i + 1) * col + j - 1] + src[(i + 1) * col + j + 1];
            // 计算垂直方向梯度
            gy = -src[(i - 1) * col + j - 1] - 2 * src[(i - 1) * col + j] - src[(i - 1) * col + j + 1] +
                 src[(i + 1) * col + j - 1] + 2 * src[(i + 1) * col + j] + src[(i + 1) * col + j + 1];
            // 根据梯度计算边缘强度
            dst[i * col + j] = sqrt(gx * gx + gy * gy) / 1141 * PIXEL_VALUE_MAX; // 1141 约为根号20，防止数据溢出
        }
    }

    return dst;
}

// 计算image中所有block的梯度直方图：基于sobel算子计算像素的梯度的大小和方向，并且根据cellSize指定block划分的cell数量
// kernel_gx = [-1 0 1 -2 0 2 -1 0 1] kernel_gy = [-1 -2 -1 0 0 0 1 2 1]
// num_bins : 角度0~180的划分bins分组数
// cell_size : 每个cell包含元素的 行/列数
// block_size : 每个block包含cell元素的 行/列数
int *calHistOfGrad(unsigned char *src, int col, int row, int num_bins, int cell_size, int block_size)
{
    if (num_bins == 0)
    {
        num_bins++; // 防止出现被除数等于0的情况
    }
    float *dst_mag = (float *)malloc(col * row * sizeof(float));
    float *dst_angle = (float *)malloc(col * row * sizeof(float));
    memset(dst_mag, 0.0, col * row * sizeof(float));
    memset(dst_angle, 0.0, col * row * sizeof(float));
    int gx, gy;
    float mag, angle;
    // 计算梯度的大小和方向角
    for (int i = 1; i < row - 1; i++)
    {
        for (int j = 1; j < col - 1; j++)
        {
            // 计算水平方向梯度
            gx = -src[(i - 1) * col + j - 1] + src[(i - 1) * col + j + 1] - 2 * src[i * col + j - 1] +
                 2 * src[i * col + j + 1] - src[(i + 1) * col + j - 1] + src[(i + 1) * col + j + 1];
            // 计算垂直方向梯度
            gy = -src[(i - 1) * col + j - 1] - 2 * src[(i - 1) * col + j] - src[(i - 1) * col + j + 1] +
                 src[(i + 1) * col + j - 1] + 2 * src[(i + 1) * col + j] + src[(i + 1) * col + j + 1];
            // 计算梯度幅值大小
            mag = sqrt(gx * gx + gy * gy) / 1141 * PIXEL_VALUE_MAX; // 1141 约为根号20，防止数据溢出
            // 计算梯度方向角
            angle = atan2(gy, gx);
            if (angle < 0)
            {
                angle += PI;
            }
            dst_mag[i * col + j] = mag;
            dst_angle[i * col + j] = angle;
        }
    }
    // 逐个Block计算梯度直方图,逐个cell滑动
    int block_row_nums = row - block_size * cell_size + 1; // 滑动block数,例 纵向共4个cell,每个block包含2*2个cell,故纵向共3个滑动block
    int block_col_nums = col - block_size * cell_size + 1;
    int *hist = (int *)malloc(block_row_nums * block_col_nums * block_size * block_size * num_bins * sizeof(int));
    memset(hist, 0, block_row_nums * block_col_nums * block_size * block_size * num_bins * sizeof(int));
    int bin, idx;
    for (int i = 0; i < block_row_nums; ++i)
    {
        for (int j = 0; j < block_col_nums; ++j)
        {
            // 统计当前Block内所有Cell中的梯度直方图
            for (int k = 0; k < block_size; ++k)
            {
                for (int t = 0; t < block_size; ++t)
                {
                    for (int m = 0; m < cell_size; ++m)
                    {
                        for (int n = 0; n < cell_size; ++n)
                        {
                            mag = dst_mag[((m + k * block_size + i) * col + n + t * block_size + j)];
                            angle = dst_angle[((m + k * block_size + i) * col + n + t * block_size + j)];
                            bin = (int)round(num_bins * angle / PI);
                            idx = (i * block_col_nums + j) * block_size * block_size * num_bins + (k * cell_size + t) * num_bins + bin;
                            hist[idx] += mag;
                        }
                    }
                }
            }
            // 对当前Block内所有Cell中的梯度直方图做向量归一化
            float block_sum = 0.0;
            for (int k = 0; k < block_size; ++k)
            {
                for (int t = 0; t < block_size; ++t)
                {
                    for (int m = 0; m < num_bins; ++m)
                    {
                        idx = (i * block_col_nums + j) * block_size * block_size * num_bins + (k * cell_size + t) * num_bins + m;
                        block_sum += hist[idx] * hist[idx];
                    }
                    block_sum = sqrt(block_sum);
                    for (int m = 0; m < num_bins; ++m)
                    {
                        idx = (i * block_col_nums + j) * block_size * block_size * num_bins + (k * cell_size + t) * num_bins + m;
                        int res = round(hist[idx] / (block_sum + 1e-6)); // +1e-6防止被除数为零的情况
                        hist[idx] = res;
                        // hist[idx] /= (block_sum - 1e-6); // +1e-6防止被除数为零的情况
                    }
                }
            }
        }
    }
    free(dst_mag);
    free(dst_angle);

    return hist;
}

void HOG(void)
{
    // 1.读取图像
    FILE *file = fopen("source_image.bmp", "rb"); // rb-以二进制模式读取
    if (file == NULL)
    {
        printf("[main err] source image file null\n");
        return;
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
    fclose(file);
    // 2.预处理 ： 灰度、归一化、均值或高斯滤波
    unsigned char *gray = rgbToGray(src, col, row);
    unsigned char *norm = normalizeImage(gray, col, row);
    unsigned char *filtered = meanFilter(norm, col, row, 2);
    // 3.边缘检测：基于sobel算子计算梯度
    unsigned char *edge = sobelDetectEgde(filtered, col, row);
    // 4.计算HOG特征

    // 内存释放
}

// int main()
// {
//     // HOG();

//     unsigned char src_arr[16] = {0, 85, 170, 255, 255, 85, 85, 0, 0, 141, 113, 0, 255, 255, 255, 255};
//     unsigned char *src = src_arr;
//     int *hist = calHistOfGrad(src, 4, 4, 9, 2, 2);

//     system("pause");

//     return 0;
// }