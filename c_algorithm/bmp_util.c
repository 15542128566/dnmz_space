#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "sift.h"

#pragma pack(1) // 结构体内存对齐设置为1字节

typedef struct BMPFILEHEADER
{
    uint16_t bfType;      // 文件类型，必须为"BM"（0x4D42）
    uint32_t bfSize;      // 文件大小，包含文件头、位图信息头和像素数据，单位为字节
    uint16_t bfReserved1; // 保留，必须为0
    uint16_t bfReserved2; // 保留，必须为0
    uint32_t bfOffBits;   // 从文件头到像素数据的偏移量，单位为字节
} BMPFILEHEADER;

typedef struct BMPINFOHEADER
{
    uint32_t biSize;         // 结构体大小，固定为40字节
    int32_t biWidth;         // 位图宽度，单位为像素
    int32_t biHeight;        // 位图高度，单位为像素，正数表示下至上，负数表示上至下
    uint16_t biPlanes;       // 颜色平面数，必须为1
    uint16_t biBitCount;     // 每像素位数，常见的为24位即RGB三色，还有32位即RGBA四色等
    uint32_t biCompression;  // 压缩方式，常见的值为0表示不压缩
    uint32_t biSizeImage;    // 像素数据大小，单位为字节，若未压缩，则值为biWidth * biHeight * biBitCount / 8
    int32_t biXPelsPerMeter; // 水平分辨率，单位为每米像素数，可设置为0表示未知
    int32_t biYPelsPerMeter; // 垂直分辨率，单位为每米像素数，可设置为0表示未知
    uint32_t biClrUsed;      // 实际使用的彩色表中的颜色索引数，为0表示默认使用所有可能的颜色
    uint32_t biClrImportant; // 对图像显示有重要影响的颜色索引数，为0表示都重要
} BMPINFOHEADER;

void save_bmp_rgb(const char *filename, uint8_t **pixels, int width, int height)
{
    BMPFILEHEADER file_header = {0x4D42, sizeof(BMPFILEHEADER) + sizeof(BMPINFOHEADER) + height * width * 3, 0, 0, sizeof(BMPFILEHEADER) + sizeof(BMPINFOHEADER)};
    BMPINFOHEADER info_header = {sizeof(BMPINFOHEADER), width, -height, 1, 24, 0, 0, 0, 0, 0, 0};

    FILE *fp = fopen(filename, "wb");
    if (!fp)
    {
        printf("Failed to open file for writing!\n");
        return;
    }

    fwrite(&file_header, 1, sizeof(file_header), fp);
    fwrite(&info_header, 1, sizeof(info_header), fp);

    int i, j;
    for (i = height - 1; i >= 0; i--)
    { // 从最后一行开始写入，遵循BMP文件格式
        for (j = 0; j < width; j++)
        {
            fwrite(pixels[i] + j * 3, 1, 3, fp); // 每个像素占用3字节，依次写入蓝色、绿色、红色
        }
    }

    fclose(fp);
}

void save_bmp_gray(int idx, unsigned char **pixels, int width, int height)
{
    BMPFILEHEADER file_header = {0x4D42, sizeof(BMPFILEHEADER) + sizeof(BMPINFOHEADER) + height * width * 3, 0, 0, sizeof(BMPFILEHEADER) + sizeof(BMPINFOHEADER)};
    // BMPFILEHEADER file_header = {0x4D42, sizeof(BMPFILEHEADER) + sizeof(BMPINFOHEADER) + height * width, 0, 0, sizeof(BMPFILEHEADER) + sizeof(BMPINFOHEADER)};
    BMPINFOHEADER info_header = {sizeof(BMPINFOHEADER), width, -height, 1, 24, 0, 0, 0, 0, 0, 0};

    char filepath[60];
    char num[5];
    char suffix[5];
    strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\siftpics\\output");
    itoa(idx, num, 10);
    strcat(filepath, num);
    strcpy(suffix, ".bmp");
    strcat(filepath, suffix);

    FILE *fp = fopen(filepath, "wb");
    if (!fp)
    {
        printf("Failed to open file for writing!\n");
        return;
    }

    fwrite(&file_header, 1, sizeof(file_header), fp);
    fwrite(&info_header, 1, sizeof(info_header), fp);

    int i, j;
    for (i = height - 1; i >= 0; i--)
    { // 从最后一行开始写入，遵循BMP文件格式
        for (j = 0; j < width; j++)
        {
            fwrite(&pixels[i][j], 1, 1, fp); // 每个像素占用3字节，依次写入蓝色、绿色、红色
            fwrite(&pixels[i][j], 1, 1, fp); // 每个像素占用3字节，依次写入蓝色、绿色、红色
            fwrite(&pixels[i][j], 1, 1, fp); // 每个像素占用3字节，依次写入蓝色、绿色、红色
            // fputc(pixels[i][j], fp); // 每个像素占用3字节，依次写入蓝色、绿色、红色
        }
    }

    fclose(fp);
}

// int main()
// {
//     int width = 320;
//     int height = 240;
//     uint8_t **pixels = (uint8_t **)malloc(height * sizeof(uint8_t *));
//     int i, j;
//     for (i = 0; i < height; i++)
//     {
//         pixels[i] = (uint8_t *)malloc(width * 3 * sizeof(uint8_t));
//         for (j = 0; j < width; j++)
//         {
//             pixels[i][j * 3 + 0] = i * 256 / height;      // 蓝色通道使用横向渐变色
//             pixels[i][j * 3 + 1] = j * 256 / width;       // 绿色通道使用纵向渐变色
//             pixels[i][j * 3 + 2] = 255 - j * 256 / width; // 红色通道使用反向渐变色
//         }
//     }
//     save_bmp("test.bmp", pixels, width, height);
//     for (i = 0; i < height; i++)
//     {
//         free(pixels[i]);
//     }
//     free(pixels);
//     return 0;
// }

void save_img_txt(unsigned char *img, int length, int idx, int type)
{
    FILE *fp;

    char filepath[60];
    char num[5];
    char suffix[5];
    if (type == 0)
    {
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_img_res\\img");
    }
    else
    {
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_pic_res\\output");
    }
    itoa(idx, num, 10);
    strcat(filepath, num);
    strcpy(suffix, ".txt");
    strcat(filepath, suffix);
    fp = fopen(filepath, "w");
    if (fp == NULL)
    {
        printf("File open error!");
        return;
    }

    for (int i = 0; i < length; i++)
    {
        fprintf(fp, "%d ", img[i]); // 将数组元素以空格分隔写入文件
    }

    fclose(fp); // 关闭文件
}

void save_mat_txt(Mat *mat, int idx, int type)
{
    FILE *fp;

    char filepath[60];
    char num[5];
    char suffix[5];
    if (type == 0)
    {
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_sample_res\\sample");
    }
    else
    {
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_pic_res\\pic");
    }
    itoa(idx, num, 10);
    strcat(filepath, num);
    strcpy(suffix, ".txt");
    strcat(filepath, suffix);
    fp = fopen(filepath, "w");
    if (fp == NULL)
    {
        printf("File open error!");
        return;
    }

    for (int i = 0; i < mat->rows; i++)
    {
        for (int j = 0; j < mat->cols; j++)
        {
            fprintf(fp, "%d ", mat->data[i][j]); // 将数组元素以空格分隔写入文件
        }
        fprintf(fp, "\n"); // 换行
    }

    fclose(fp); // 关闭文件
}

void save_double_mat_txt(DoubleMat *double_mat, int idx, int type)
{
    FILE *fp;

    char filepath[60];
    char num[5];
    char suffix[5];
    switch (type)
    {
    case 0:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_guas_res\\guas");
        break;
    case 1:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_kernel_res\\kernel");
        break;
    case 2:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_diff_res\\diff");
        break;
    default:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_guas_res\\guas");
        break;
    }
    itoa(idx, num, 10);
    strcat(filepath, num);
    strcpy(suffix, ".txt");
    strcat(filepath, suffix);
    fp = fopen(filepath, "w");
    if (fp == NULL)
    {
        printf("File open error!");
        return;
    }

    for (int i = 0; i < double_mat->rows; i++)
    {
        for (int j = 0; j < double_mat->cols; j++)
        {
            fprintf(fp, "%f ", double_mat->data[i][j]); // 将数组元素以空格分隔写入文件
        }
        fprintf(fp, "\n"); // 换行
    }

    fclose(fp); // 关闭文件
}

void save_key_points_txt(PointList *point_list, int type)
{
    FILE *fp;

    char filepath[60];
    char num[5];
    char suffix[5];
    switch (type)
    {
    case 0:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\key_points");
        break;
    default:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\key_points");
        break;
    }
    strcpy(suffix, ".txt");
    strcat(filepath, suffix);
    fp = fopen(filepath, "w");
    if (fp == NULL)
    {
        printf("File open error!");
        return;
    }
    Point point;
    for (int i = 0; i < point_list->length; i++)
    {
        point = point_list->point_list[i];
        fprintf(fp, "%f ", point.posX); // 将数组元素以空格分隔写入文件
        fprintf(fp, "%f ", point.posY);
        fprintf(fp, "%d ", point.posOS);
        fprintf(fp, "%f ", point.imgScale);
        fprintf(fp, "%f ", point.direction);
        fprintf(fp, "\n"); // 换行
    }

    fclose(fp); // 关闭文件
}

void save_descriptor_txt(DescriptorList *dpr_list, int type)
{
    FILE *fp;

    char filepath[60];
    char num[5];
    char suffix[5];
    switch (type)
    {
    case 0:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\descriptor");
        break;
    default:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\descriptor");
        break;
    }
    strcpy(suffix, ".txt");
    strcat(filepath, suffix);
    fp = fopen(filepath, "w");
    if (fp == NULL)
    {
        printf("File open error!");
        return;
    }
    Descriptor dpr;
    for (int i = 0; i < dpr_list->length; i++)
    {
        dpr = dpr_list->descriptor_list[i];
        for (int j = 0; j < dpr.length; j++)
        {
            fprintf(fp, "%f ", dpr.data[j]); // 将数组元素以空格分隔写入文件
        }
        fprintf(fp, "\n"); // 换行
    }

    fclose(fp); // 关闭文件
}

void save_descriptor_tmp_txt(double *data, int len, int type)
{
    FILE *fp;

    char filepath[60];
    char num[5];
    char suffix[5];
    switch (type)
    {
    case 0:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\x_list");
        break;
    case 1:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\y_list");
        break;
    case 2:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\rbin_list");
        break;
    case 3:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\cbin_list");
        break;
    case 4:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\w_list");
        break;
    case 5:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\mag_ist");
        break;
    case 6:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\angle_list");
        break;
    case 7:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\hist_list");
        break;
    default:
        strcpy(filepath, "C:\\Users\\A\\Desktop\\sift\\sift_descriptor_res\\descriptor");
        break;
    }
    strcpy(suffix, ".txt");
    strcat(filepath, suffix);
    fp = fopen(filepath, "w");
    if (fp == NULL)
    {
        printf("File open error!");
        return;
    }
    Descriptor dpr;
    for (int i = 0; i < len; i++)
    {
        fprintf(fp, "%f ", data[i]); // 将数组元素以空格分隔写入文件
    }
    fprintf(fp, "\n"); // 换行

    fclose(fp); // 关闭文件
}