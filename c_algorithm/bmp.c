#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#pragma pack(2)   //以 2 个字节对齐
typedef struct tagBITMAPFILEHEADER{
    unsigned short    bfType;
    unsigned long     bfSize;
    unsigned short    bfReserved1;
    unsigned short    bfReserved2;
    unsigned long     bfOffBits;
}BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER{
    unsigned long     biSize;
    long              biWidth;
    long              biHeight;
    unsigned short    biPlanes;
    unsigned short    biBitCount;
    unsigned long     biCompression;
    unsigned long     biSizeImage;
    long              biXPelsPerMeter;
    long              biYPelsPerMeter;
    unsigned long     biClrUsed;
    unsigned long     biClrImportant;
}BITMAPINFOHEADER;

#pragma pack()

// 从文件中读取一个 unsigned char 类型的数据
unsigned char read_uchar(FILE *fp)
{
    unsigned char b;
    fread(&b, 1, 1, fp);
    return b;
}

// 从文件中读取一个 unsigned short 类型的数据
unsigned short read_ushort(FILE *fp)
{
    unsigned short s;
    fread(&s, 2, 1, fp);
    return s;
}

// 从文件中读取一个 unsigned int 类型的数据
unsigned int read_uint(FILE *fp)
{
    unsigned int i;
    fread(&i, 4, 1, fp);
    return i;
}

// 从文件中读取一个 int 类型的数据
int read_int(FILE *fp)
{
    int i;
    fread(&i, 4, 1, fp);
    return i;
}

// 从文件中读取 BITMAPFILEHEADER 结构体数据
void read_bitmap_file_header(FILE *fp, BITMAPFILEHEADER *pFileHeader)
{
    pFileHeader->bfType          = read_ushort(fp);
    pFileHeader->bfSize          = read_uint(fp);
    pFileHeader->bfReserved1     = read_ushort(fp);
    pFileHeader->bfReserved2     = read_ushort(fp);
    pFileHeader->bfOffBits       = read_uint(fp);
}

// 从文件中读取 BITMAPINFOHEADER 结构体数据
void read_bitmap_info_header(FILE *fp, BITMAPINFOHEADER *pInfoHeader)
{
    pInfoHeader->biSize          = read_uint(fp);
    pInfoHeader->biWidth         = read_int(fp);
    pInfoHeader->biHeight        = read_int(fp);
    pInfoHeader->biPlanes        = read_ushort(fp);
    pInfoHeader->biBitCount      = read_ushort(fp);
    pInfoHeader->biCompression   = read_uint(fp);
    pInfoHeader->biSizeImage     = read_uint(fp);
    pInfoHeader->biXPelsPerMeter = read_int(fp);
    pInfoHeader->biYPelsPerMeter = read_int(fp);
    pInfoHeader->biClrUsed       = read_uint(fp);
    pInfoHeader->biClrImportant  = read_uint(fp);
}

int main()
{
    // 打开 BMP 文件
    FILE* fp = fopen("test.bmp", "rb");
    if(fp == NULL){
        printf("File open failed!\n");
        return -1;
    }

    BITMAPFILEHEADER fileHeader;
    BITMAPINFOHEADER infoHeader;

    // 读取 BMP 文件头数据
    read_bitmap_file_header(fp, &fileHeader);

    // 读取 BMP 图像信息头数据，获取图像的宽度、高度、位数等信息
    read_bitmap_info_header(fp, &infoHeader);
    printf("width: %d, height: %d, bitCount: %d\n",
        infoHeader.biWidth, infoHeader.biHeight, infoHeader.biBitCount);

    // 计算 BMP 图像数据区域的大小
    unsigned int dataSize = infoHeader.biSizeImage != 0 ? infoHeader.biSizeImage :
                            infoHeader.biWidth * infoHeader.biHeight * (infoHeader.biBitCount / 8);

    // 根据图像数据区域的大小，动态分配一个缓冲区
    unsigned char* buffer = (unsigned char*)malloc(dataSize);
    if(buffer == NULL){
        fclose(fp);
        return -1;
    }

    // 将指针移动到图像数据的起始位置
    fseek(fp, fileHeader.bfOffBits, SEEK_SET);

    // 读取图像数据到缓冲区
    fread(buffer, 1, dataSize, fp);

    fclose(fp);

    // 从缓冲区中读取像素点数据
    int bytePerPix = infoHeader.biBitCount / 8;   //每个像素有多少字节
    int width = infoHeader.biWidth;
    int height = abs(infoHeader.biHeight);
    unsigned char* pData = buffer;
    for(int y = 0; y < height; y++){
        for(int x = 0; x < width; x++){
            unsigned char b = *pData;
            unsigned char g = *(pData + 1);
            unsigned char r = *(pData + 2);
            printf("x: %d, y: %d, r: %d, g: %d, b: %d\n", x, y, r, g, b);
            pData += bytePerPix;
        }
    }

    free(buffer);

    return 0;
}
