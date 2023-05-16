#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RGB_COUNT 3         // 每个像素包含色彩数
#define PIXEL_VALUE_MIN 0   // 像素最小范围
#define PIXEL_VALUE_MAX 255 // 像素最大范围

static const int BMP_FILE_HEADER_LEN = 14;                                     // 位图文件BMP格式，文件头14字节
static const int BMP_INFO_HEADER_LEN = 40;                                     // 位图文件BMP格式，信息头40字节
static const float COF_R2GRAY = 0.299, COF_G2GRAY = 0.587, COF_B2GRAY = 0.114; // rgb转化为灰度像素的系数

// rgbToGray将图像处理为灰度图像
// 共row*col个像素，每个像素包含rgb3个字节，逐行逐列个像素处理
unsigned char *rgbToGray(unsigned char *src, int row, int col)
{
    unsigned char *gray = (unsigned char *)malloc(row * col * sizeof(unsigned char));
    for (int i = 0; i < col; i++)
    {
        for (int j = 0; j < row; j++)
        {
            int idx = (i * row + j) * RGB_COUNT;
            int r = src[idx];
            int g = src[idx + 1];
            int b = src[idx + 2];
            gray[i * row + j] = COF_R2GRAY * r + COF_G2GRAY * g + COF_B2GRAY * b;
        }
    }

    return gray;
}

// normalizeImage归一化处理，将图像中的每个像素缩放至[0~255]范围
unsigned char *normalizeImage(unsigned char *src, int row, int col)
{
    unsigned char *dst = (unsigned char *)malloc(row * col * sizeof(unsigned char));
    // 找到所有像素点中的最大最小值, 计算缩放比例scale
    int min = PIXEL_VALUE_MAX;
    int max = PIXEL_VALUE_MIN;
    for (int i = 0; i < col; i++)
    {
        for (int j = 0; j < row; j++)
        {
            int pv = src[i * row + j];
            min = (pv < min) ? pv : min;
            max = (pv > max) ? pv : max;
        }
    }
    if (max > min)
    {
        double scale = PIXEL_VALUE_MAX / (max - min);
        for (int i = 0; i < col; i++)
        {
            for (int j = 0; j < row; j++)
            {
                src[i * row + j] = (unsigned char)((src[i * row + j] - min) * scale);
            }
        }
    }
    else
    {
        printf("[normalizeImage warn] max value < min value");
        memcpy(dst, src, row * col); // 若max<min则放弃归一化，直接返回原始值
    }

    return dst;
}

// 均值滤波(边缘做截断处理) winsize为窗口边数
unsigned char *averageFilter(unsigned char *src, int row, int col, int winsize)
{
    if (row < winsize || col < winsize)
    {
        printf("averageFilter error: invalid params, row:%d,col:%d,winsize:%d\n", row, col, winsize);

        return NULL;
    }
    if (winsize & 1)
    {
        winsize += 1; //窗口边数调整为奇数
    }
    int winsize_2 = winsize/2;
    unsigned char *dst = (unsigned char *)malloc(row * col * sizeof(unsigned char));
    double average;
    int idx;
    for (int i = winsize_2; i < col - winsize_2; i++)
    {
        for (int j = winsize_2; j < row - winsize_2; j++)
        {
            average = 0.0;
            for (int k = 0; k < winsize; k++)
            {
                for (int s = 0; s < winsize; s++)
                {
                    idx =
                }
            }
        }
    }
}

void HOG(void)
{
    // 1.读取图像
    FILE *file = fopen("source_image.bmp", "rb"); // rb-以二进制模式读取
    if (file == NULL)
    {
        printf("[main err] source image file null");
        return 0;
    }
    // 获取图像大小
    int row = 0, col = 0;                           // 图像宽高
    fseek(file, BMP_FILE_HEADER_LEN + 4, SEEK_SET); // 跳过文件头和信息头的第一个字段
    fread(&row, sizeof(int), 1, file);
    fread(&col, sizeof(int), 1, file);
    // 存储图像数据
    unsigned char *src = NULL;
    src = (unsigned char *)malloc(row * col * RGB_COUNT * sizeof(unsigned char));
    fseek(file, (BMP_FILE_HEADER_LEN + BMP_INFO_HEADER_LEN), SEEK_SET); // 跳过文件头和信息头
    // todo: 应该根据颜色表判断的存储方式为RGB还是索引方式，这里简化处理为默认RGB存储方式
    fread(src, sizeof(char), row * col * RGB_COUNT, file);
    fclose(file);
    // 2.预处理 ： 灰度、归一化、均值或高斯滤波、边缘检测

    // 3.计算HOG特征

    // 内存释放
}

int main()
{
    HOG();

    return 0;
}