#include <gtest/gtest.h>
#include "../c_algorithm/hog.c"

// 测试gtest
TEST(ExampleTests, DemonstrateGTestMacros)
{
    EXPECT_TRUE(true) << "Hello, world!";
    EXPECT_EQ(true, true) << "Fuck, world!";
}

// 测试图像像素灰度化函数 rgbToGray
TEST(testGray, rgbToGray)
{
    unsigned char rgb[16 * 3] = {1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
                                 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1,
                                 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3,
                                 4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6};
    unsigned char *src = rgb;
    unsigned char *gray = rgbToGray(src, 4, 4);
    unsigned char expectexpectGrayArr[16] = {1, 2, 3, 4, 4, 3, 2, 1, 1, 1, 1, 1, 4, 4, 4, 4};
    unsigned char *expectGray = expectexpectGrayArr;
    for (int i = 0; i < 16; ++i)
    {
        EXPECT_EQ(
            expectGray[i],
            gray[i])
            << "i=" << i << ",gray_i=" << static_cast<int>(gray[i]);
    };
}

// 测试图像像素归一化函数 normalizeImage
TEST(testNormalize, normalizeImage)
{
    unsigned char src_arr[16] = {1, 2, 3, 4, 4, 3, 2, 1, 1, 1, 1, 1, 4, 4, 4, 4};
    unsigned char *src = src_arr;
    unsigned char *norm = normalizeImage(src, 4, 4);
    unsigned char expertNormArr[16] = {0, 85, 170, 255, 255, 170, 85, 0, 0, 0, 0, 0, 255, 255, 255, 255};
    unsigned char *expectNorm = expertNormArr;
    for (int i = 0; i < 16; ++i)
    {
        EXPECT_EQ(
            expectNorm[i],
            norm[i])
            << "i=" << i << ",norm_i=" << static_cast<int>(norm[i]);
    };
}

// 测试均值滤波处理 meanFilter
TEST(testMeanFilter, meanFilter)
{
    unsigned char src_arr[16] = {0, 85, 170, 255, 255, 170, 85, 0, 0, 0, 0, 0, 255, 255, 255, 255};
    unsigned char *src = src_arr;
    unsigned char *mean = meanFilter(src, 4, 4, 2);
    unsigned char expectMeanArr[16] = {0, 85, 170, 255, 255, 85, 85, 0, 0, 141, 113, 0, 255, 255, 255, 255};
    unsigned char *expectMean = expectMeanArr;
    for (int i = 0; i < 16; ++i)
    {
        EXPECT_EQ(
            expectMean[i],
            mean[i])
            << "i=" << i << ",mean_i=" << static_cast<int>(mean[i]);
    };
}

// 测试基于sobel梯度算子的图像边缘检测
TEST(testSobelDetectEdge, sobelDetectEgde)
{
    unsigned char src_arr[16] = {0, 85, 170, 255, 255, 85, 85, 0, 0, 141, 113, 0, 255, 255, 255, 255};
    unsigned char *src = src_arr;
    unsigned char *sobel = sobelDetectEgde(src, 4, 4);
    unsigned char expectSobelArr[16] = {0, 0, 0, 0, 0, 17, 76, 0, 0, 114, 189, 0, 0, 0, 0, 0};
    unsigned char *expectSobel = expectSobelArr;
    for (int i = 0; i < 16; ++i)
    {
        EXPECT_EQ(
            expectSobel[i],
            sobel[i])
            << "i=" << i
            << ",sobel_i=" << static_cast<int>(sobel[i]);
    };
}

// 测试基于sobel梯度算子的梯度直方图特征向量
TEST(testcalHistOfGrad, calHistOfGrad)
{
    unsigned char src_arr[16] = {0, 85, 170, 255, 255, 85, 85, 0, 0, 141, 113, 0, 255, 255, 255, 255};
    unsigned char *src = src_arr;
    int *hist = calHistOfGrad(src, 4, 4, 9, 2, 2);
    int expectHistArr[36] = {0, 0, 0, 0, 0, 0, 0, 1, 0,
                             0, 0, 0, 1, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 1, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 1, 0, 0};
    int *expectHist = expectHistArr;
    for (int i = 0; i < 36; ++i)
    {
        EXPECT_EQ(
            expectHist[i],
            hist[i])
            << "i=" << i
            << ",hist_i=" << hist[i];
    };
}